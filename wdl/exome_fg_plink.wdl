version 1.0

# Subset FinnGen VCFs to exome samples and convert to plink per chromosome.
#
# Workflow shape:
#   QueryChromPositions       single task — all 23 chroms queried in parallel (cached)
#   BuildFGChunkRegions       single task — extracts sample lists on the fly, then
#                             builds chunk regions bounded by output size and
#                             input streaming time (token safety)
#   SubsetFGChunkToPlink[K]   one task per chunk — bcftools view → tmp.bcf → plink2
#   GatherChromChunks         single task — groups flat chunk paths by chrom (no
#                             localization: Array[String] in, Array[Array[File]] out)
#   MergeFGChromChunks[C]     one task per chrom — plink2 --pmerge-list

workflow exome_fg_plink {
  input {
    Array[Array[String]] exome_vcf_pairs   # [[prefix, gs://…_chrCHROM.vcf.gz], …]
    String               fg_vcf_basepath   # gs://…_chrCHROM.vcf.gz  (CHROM = 1, X, …)
    Array[String]        chroms            # e.g. ["chr1", "chrX"]
    Int                  chunk_mb  = 500   # target output plink chunk size in MB
    String               plink_conv_args   # e.g. "--double-id --allow-extra-chr …"
    String               plink_merge_args  # e.g. "--allow-extra-chr"
  }

  # ── PHASE 1: query FG positions for all chroms in parallel (cached) ─────────
  call QueryChromPositions {
    input:
      fg_vcf_basepath = fg_vcf_basepath
  }

  # ── PHASE 2: extract samples on the fly + build chunk regions ──────────────
  call BuildFGChunkRegions {
    input:
      all_positions   = QueryChromPositions.positions,
      exome_vcf_pairs = exome_vcf_pairs,
      fg_vcf_basepath = fg_vcf_basepath,
      chroms          = chroms,
      chunk_mb        = chunk_mb
  }

  # ── PHASE 3: subset + convert each chunk to plink ──────────────────────────
  scatter (i in range(length(BuildFGChunkRegions.tasks))) {
    call SubsetFGChunkToPlink {
      input:
        out_prefix      = BuildFGChunkRegions.tasks[i][0],
        fg_vcf          = BuildFGChunkRegions.tasks[i][1],
        region          = BuildFGChunkRegions.tasks[i][2],
        union_samples   = BuildFGChunkRegions.union_samples,
        plink_conv_args = plink_conv_args,
        chunk_mb        = chunk_mb
    }
  }

  # ── PHASE 4: group flat chunk paths by chrom (no localization) ────────────
  call GatherChromChunks {
    input:
      all_beds = SubsetFGChunkToPlink.out_bed,
      all_bims = SubsetFGChunkToPlink.out_bim,
      all_fams = SubsetFGChunkToPlink.out_fam,
      chroms   = chroms
  }

  # ── PHASE 5: merge chunks → one plink per chrom ────────────────────────────
  scatter (ch_i in range(length(chroms))) {
    call MergeFGChromChunks {
      input:
        chrom_files      = GatherChromChunks.chrom_files[ch_i],
        chrom            = chroms[ch_i],
        chunk_mb         = chunk_mb,
        plink_merge_args = plink_merge_args
    }
  }

  output {
    Array[File] fg_beds = MergeFGChromChunks.out_bed
    Array[File] fg_bims = MergeFGChromChunks.out_bim
    Array[File] fg_fams = MergeFGChromChunks.out_fam
  }
}


# ---------------------------------------------------------------------------
# Query POS for all 23 chroms in parallel. Each command refreshes its own
# token so long-running queries don't expire. Runs once and is cached.
# ---------------------------------------------------------------------------
task QueryChromPositions {
  input {
    String fg_vcf_basepath
    Int    cpu     = 23
    Int    disk_gb = 30
  }

  command <<<
  set -euo pipefail

  for chrom in {1..22} X; do
    vcf=$(echo "~{fg_vcf_basepath}" | sed "s/CHROM/$chrom/")
    echo "export GCS_OAUTH_TOKEN=\$(gcloud auth application-default print-access-token) && \
      bcftools query -f '%POS\n' \"$vcf\" | gzip > chr${chrom}_pos.txt.gz && \
      echo \"chr${chrom}: \$(zcat chr${chrom}_pos.txt.gz | wc -l) positions\" >&2"
  done > query_commands.sh

  parallel -j ~{cpu} < query_commands.sh
  >>>

  output {
    Array[File] positions = glob("chr*_pos.txt.gz")
  }

  runtime {
    cpu:   cpu
    disks: "local-disk ~{disk_gb} HDD"
  }
}


# ---------------------------------------------------------------------------
# Extract sample lists on the fly (header-only, seconds per VCF), then build
# size-based chunk regions for all chroms.
#
# Chunk size is the minimum of two bounds:
#   output bound  — expected plink output after sample subsetting <= chunk_mb
#   input bound   — input bytes to stream from GCS <= ~10 GB (token safety)
#
# Output tasks TSV columns: out_prefix | vcf_path | region
# ---------------------------------------------------------------------------
task BuildFGChunkRegions {
  input {
    Array[File]          all_positions
    Array[Array[String]] exome_vcf_pairs
    String               fg_vcf_basepath
    Array[String]        chroms
    Int                  chunk_mb
    Int                  disk_gb = 20
  }

  command <<<
  set -euo pipefail
  export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

  # Extract and merge exome sample lists (header-only reads, seconds each)
  while IFS=$'\t' read -r prefix vcf_base; do
    chr1_vcf=$(echo "$vcf_base" | sed 's/CHROM/1/')
    bcftools query -l "$chr1_vcf"
  done < ~{write_tsv(exome_vcf_pairs)} | sort -u > union_samples.txt
  n_exome=$(wc -l < union_samples.txt | tr -d ' ')

  # Get FG sample count from chr1 header
  fg_chr1=$(echo "~{fg_vcf_basepath}" | sed 's/CHROM/1/')
  n_fg=$(bcftools query -l "$fg_chr1" | wc -l | tr -d ' ')

  echo "exome union: $n_exome  FG: $n_fg" >&2

  # Match each requested chrom to its position file by basename
  while read chrom; do
    pos_file=$(grep "/${chrom}_pos\.txt\.gz" ~{write_lines(all_positions)})
    printf "%s\t%s\n" "$chrom" "$pos_file"
  done < ~{write_lines(chroms)} > chrom_pos_pairs.txt

  while IFS=$'\t' read -r chrom pos_file; do
    chrom_num=$(echo "$chrom" | sed 's/chr//')
    vcf_path=$(echo "~{fg_vcf_basepath}" | sed "s/CHROM/$chrom_num/")

    file_bytes=$(gsutil du "$vcf_path" | awk '{print $1}')
    n_variants=$(bcftools index --stats "$vcf_path" | awk -v c="$chrom" '$1==c{print $3}')

    input_bpv=$(awk  -v fb="$file_bytes" -v nv="$n_variants" \
      'BEGIN{printf "%.4f", fb/nv}')
    output_bpv=$(awk -v ibpv="$input_bpv" -v ne="$n_exome" -v nf="$n_fg" \
      'BEGIN{printf "%.4f", ibpv*(ne/nf)}')

    n_per_chunk=$(awk \
      -v obpv="$output_bpv" -v ibpv="$input_bpv" -v mb="~{chunk_mb}" \
      'BEGIN{
        n_out = int(mb * 1024 * 1024 / obpv);
        n_in  = int(10  * 1024 * 1024 * 1024 / ibpv);
        n     = (n_out < n_in ? n_out : n_in);
        print (n < 1 ? 1 : n)
      }')

    echo "$chrom: input_bpv=$input_bpv output_bpv=$output_bpv n_per_chunk=$n_per_chunk n_variants=$n_variants" >&2

    last_pos=$(zcat "$pos_file" | tail -n1)
    { zcat "$pos_file"; echo "$last_pos"; } | split -l "$n_per_chunk" - "split_${chrom}_"

    chunk=0
    start=1
    for f in $(ls "split_${chrom}_"* | sort); do
      end=$(tail -n1 "$f")
      printf "%s_%03d\t%s\t%s:%d-%d\n" \
        "$chrom" "$chunk" "$vcf_path" "$chrom" "$start" "$end"
      start=$((end + 1))
      chunk=$((chunk + 1))
    done
    rm "split_${chrom}_"*
  done < chrom_pos_pairs.txt > tasks.tsv

  echo "Total chunks: $(wc -l < tasks.tsv)" >&2
  >>>

  output {
    Array[Array[String]] tasks         = read_tsv("tasks.tsv")
    File                 union_samples = "union_samples.txt"
  }

  runtime {
    disks: "local-disk ~{disk_gb} HDD"
  }
}


# ---------------------------------------------------------------------------
# Stream one FG region, subset to exome samples, convert to plink.
# ---------------------------------------------------------------------------
task SubsetFGChunkToPlink {
  input {
    String out_prefix
    String fg_vcf
    String region
    File   union_samples
    String plink_conv_args
    Int    chunk_mb
    Int    cpu     = 8
    Int    disk_gb = chunk_mb * 3 / 1024 + 10
  }

  command <<<
  set -euo pipefail
  export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
  export HTS_HTTP_VERSION=1.1

  bcftools view \
    --regions       "~{region}" \
    --samples-file  "~{union_samples}" \
    --force-samples \
    --threads       ~{cpu} \
    --output-type   b \
    --output        tmp.bcf \
    "~{fg_vcf}"

  plink2 \
    --bcf       tmp.bcf \
    --make-bed \
    --out       "~{out_prefix}" \
    --threads   ~{cpu} \
    ~{plink_conv_args}

  rm tmp.bcf

  echo "DONE ~{out_prefix}" \
    "variants=$(wc -l < ~{out_prefix}.bim)" \
    "samples=$(wc -l < ~{out_prefix}.fam)" >&2
  >>>

  output {
    File out_bed = out_prefix + ".bed"
    File out_bim = out_prefix + ".bim"
    File out_fam = out_prefix + ".fam"
  }

  runtime {
    cpu:   cpu
    disks: "local-disk ~{disk_gb} HDD"
  }
}


# ---------------------------------------------------------------------------
# Group flat chunk path arrays by chrom. Inputs are String (no localization);
# bed/bim/fam are combined into a single row per chrom so the merge task
# receives one Array[File] with all three types — plink finds them by prefix.
# ---------------------------------------------------------------------------
task GatherChromChunks {
  input {
    Array[String] all_beds
    Array[String] all_bims
    Array[String] all_fams
    Array[String] chroms
    Int           disk_gb = 5
  }

  command <<<
  set -euo pipefail

  while read chrom; do
    {
      grep "${chrom}_[0-9]" ~{write_lines(all_beds)} | awk -F/ '{print $NF"\t"$0}' | sort -k1,1 | cut -f2-
      grep "${chrom}_[0-9]" ~{write_lines(all_bims)} | awk -F/ '{print $NF"\t"$0}' | sort -k1,1 | cut -f2-
      grep "${chrom}_[0-9]" ~{write_lines(all_fams)} | awk -F/ '{print $NF"\t"$0}' | sort -k1,1 | cut -f2-
    } | paste -s -d'\t' >> chrom_files.tsv
  done < ~{write_lines(chroms)}
  >>>

  output {
    Array[Array[File]] chrom_files = read_tsv("chrom_files.tsv")
  }

  runtime {
    disks: "local-disk ~{disk_gb} HDD"
  }
}


# ---------------------------------------------------------------------------
# Merge all plink chunks for one chrom via plink2 --pmerge-list.
# Receives bed+bim+fam for this chrom only — Cromwell localizes exactly those.
# ---------------------------------------------------------------------------
task MergeFGChromChunks {
  input {
    Array[File] chrom_files
    String      chrom
    Int         chunk_mb
    String      plink_merge_args
    Int         cpu     = 8
    Int         disk_gb = chunk_mb * 3 / 1024 + 10
  }

  String out_prefix = "fg_" + chrom

  command <<<
  set -euo pipefail

  for f in ~{sep=" " chrom_files}; do
    ln -sf "$f" "$(basename $f)"
  done

  ls *.bed | sed 's/\.bed$//' | sort > merge_list.txt
  n_chunks=$(wc -l < merge_list.txt)
  echo "~{chrom}: merging $n_chunks chunks" >&2

  if [[ $n_chunks -eq 1 ]]; then
    prefix=$(cat merge_list.txt)
    mv "${prefix}.bed" "~{out_prefix}.bed"
    mv "${prefix}.bim" "~{out_prefix}.bim"
    mv "${prefix}.fam" "~{out_prefix}.fam"
  else
    plink2 \
      --pmerge-list merge_list.txt bfile \
      --make-bed \
      --out       "~{out_prefix}" \
      --threads   ~{cpu} \
      ~{plink_merge_args}
  fi

  echo "DONE ~{out_prefix}" \
    "variants=$(wc -l < ~{out_prefix}.bim)" \
    "samples=$(wc -l < ~{out_prefix}.fam)" >&2
  >>>

  output {
    File out_bed = out_prefix + ".bed"
    File out_bim = out_prefix + ".bim"
    File out_fam = out_prefix + ".fam"
  }

  runtime {
    cpu:   cpu
    disks: "local-disk ~{disk_gb} HDD"
  }
}
