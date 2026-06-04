version 1.0

# Subsample and rename exome VCFs using the QRY→REF mapping produced by
# exome_duplicates.wdl.
#
# Workflow shape:
#   QueryChromPositions[C]  one task per chrom — queries positions for all
#                           datasets in parallel (cached per chrom)
#   BuildAllRegions         single task — uses position files + file stats to
#                           build per-dataset size-based chunk TSV
#   SubsetChunk[K]          one task per (dataset × chrom × chunk)
#   ConcatChromVCF[D×C]     cross-product scatter — greps flat VCF path array
#                           to concat chunks for one dataset+chrom

workflow exome_rename {
  input {
    File                 resolved_mapping
    Array[Array[String]] vcf_pairs = [
      ["BOTNIA", "gs://fg-3/exome_v2/data/THLBB2023_14_WES_Botnia.QC_ANNOTATED.vcf.gz"],
      ["ADPKD",  "gs://fg-3/exome_v2/data/likely_pathogenic_annot_annotated_full_header_fix_resampled.QC_ANNOTATED.vcf.gz"],
      ["DALY",   "gs://fg-3/exome_v2/data/fimm-daly_finnish_gvs_bge_callset_1_padded_split_FINBBonly.QC_ANNOTATED.vcf.gz"],
      ["WES",    "gs://fg-3/exome_v2/data/finngen_wes_gnomad_v4.QC_ANNOTATED.vcf.gz"]
    ]
    Array[String] chroms = [
      "chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
      "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
      "chr20","chr21","chr22","chrX","chrY"
    ]
    Int    chunk_mb = 500
    String suffix   = "fg_ids"
  }

  scatter (chrom in chroms) {
    call QueryChromPositions {
      input:
        chrom     = chrom,
        vcf_pairs = vcf_pairs
    }
  }

  call BuildAllRegions {
    input:
      position_files = flatten(QueryChromPositions.positions),
      vcf_pairs      = vcf_pairs,
      chroms         = chroms,
      chunk_mb       = chunk_mb
  }

  scatter (i in range(length(BuildAllRegions.tasks))) {
    call SubsetChunk {
      input:
        vcf              = BuildAllRegions.tasks[i][0],
        prefix           = BuildAllRegions.tasks[i][1],
        region           = BuildAllRegions.tasks[i][2],
        resolved_mapping = resolved_mapping,
        chunk_mb         = chunk_mb
    }
  }

  Int n_gather = length(vcf_pairs) * length(chroms)

  scatter (i in range(n_gather)) {
    Int ds_idx = i / length(chroms)
    Int ch_idx = i % length(chroms)
    call ConcatChromVCF {
      input:
        all_vcf_paths = SubsetChunk.out_vcf,
        dataset       = vcf_pairs[ds_idx][0],
        vcf           = vcf_pairs[ds_idx][1],
        chrom         = chroms[ch_idx],
        suffix        = suffix,
        chrom_size_mb = BuildAllRegions.chrom_sizes_mb[i]
    }
  }

  output {
    Array[File] chrom_vcfs = ConcatChromVCF.out_vcf
    Array[File] chrom_tbis = ConcatChromVCF.out_tbi
  }
}


# ---------------------------------------------------------------------------
# Query CHROM+POS for all datasets for one chrom in parallel (cached).
# Output: one {dataset}_{chrom}_pos.txt per dataset.
# ---------------------------------------------------------------------------
task QueryChromPositions {
  input {
    String               chrom
    Array[Array[String]] vcf_pairs
    Int                  cpu     = 4
    Int                  disk_gb = 10
  }

  command <<<
  set -euo pipefail

  while IFS=$'\t' read -r dataset vcf; do
    echo "export GCS_OAUTH_TOKEN=\$(gcloud auth application-default print-access-token) && bcftools query -f '%POS\n' --regions ~{chrom} \"$vcf\" > ${dataset}_~{chrom}_pos.txt"
  done < ~{write_tsv(vcf_pairs)} > query_commands.sh

  parallel -j ~{cpu} < query_commands.sh

  for f in *_~{chrom}_pos.txt; do
    echo "~{chrom} $(basename $f _~{chrom}_pos.txt): $(wc -l < $f) positions" >&2
  done
  >>>

  output {
    Array[File] positions = glob("*_~{chrom}_pos.txt")
  }

  runtime {
    cpu:   cpu
    disks: "local-disk ~{disk_gb} HDD"
  }
}


# ---------------------------------------------------------------------------
# Build size-based chunk regions for all datasets × chroms.
# Queries file stats (fast — index only), then loops over position files.
# Output: Array[Array[String]] tasks with columns vcf | prefix | region
# prefix = DATASET_CHUNKIDX e.g. BOTNIA_001
# ---------------------------------------------------------------------------
task BuildAllRegions {
  input {
    Array[File]          position_files
    Array[Array[String]] vcf_pairs
    Array[String]        chroms
    Int                  chunk_mb
  }

  command <<<
  set -euo pipefail
  export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

  # Link position files by name
  while read -r f; do
    ln -s "$f" "$(basename $f)"
  done < ~{write_lines(position_files)}

  # Get file size + total variants per dataset
  while IFS=$'\t' read -r dataset vcf; do
    file_bytes=$(gsutil du "$vcf" | awk '{print $1}')
    n_total=$(bcftools index --stats "$vcf" | awk '{sum+=$3} END{print sum}')
    bpv=$(awk -v fb="$file_bytes" -v nt="$n_total" 'BEGIN{printf "%.2f", fb/nt}')
    echo "$dataset $vcf $bpv" >> dataset_stats.txt
    echo "$dataset: $(( file_bytes/1024/1024 ))MB  total_variants=$n_total  bpv=$bpv" >&2
  done < ~{write_tsv(vcf_pairs)}

  # Build regions per dataset × chrom; also emit chrom_sizes_mb in scatter order
  while read -r dataset vcf bpv; do
    while read -r chrom; do
      pos_file="${dataset}_${chrom}_pos.txt"
      if [[ ! -f "$pos_file" || ! -s "$pos_file" ]]; then
        echo "0" >> chrom_sizes_mb.txt
        continue
      fi

      n_chrom=$(wc -l < "$pos_file")
      chrom_mb=$(awk -v bpv="$bpv" -v nc="$n_chrom" 'BEGIN{printf "%.0f", bpv*nc/1024/1024}')
      echo "$chrom_mb" >> chrom_sizes_mb.txt

      n_variants=$(awk -v bpv="$bpv" -v mb="~{chunk_mb}" 'BEGIN{n=int(mb*1024*1024/bpv); print (n<1?1:n)}')
      last_pos=$(tail -n1 "$pos_file")

      { cat "$pos_file"; echo "$last_pos"; } | split -l "$n_variants" - "split_${dataset}_${chrom}_"

      chunk=0
      start=1
      for f in $(ls split_${dataset}_${chrom}_* | sort); do
        end=$(tail -n1 "$f")
        printf "%s\t%s_%s_%03d\t%s:%d-%d\n" "$vcf" "$dataset" "$chrom" "$chunk" "$chrom" "$start" "$end"
        start=$((end + 1))
        chunk=$((chunk + 1))
      done
      rm split_${dataset}_${chrom}_*
    done < ~{write_lines(chroms)}
  done < dataset_stats.txt > tasks.tsv

  echo "Total scatter tasks: $(wc -l < tasks.tsv)" >&2
  >>>

  output {
    Array[Array[String]] tasks          = read_tsv("tasks.tsv")
    Array[Int]           chrom_sizes_mb = read_lines("chrom_sizes_mb.txt")
  }

  runtime {
    disks: "local-disk 20 HDD"
  }
}


# ---------------------------------------------------------------------------
# Stream one region, subset + rename samples.
# prefix = DATASET_CHUNKIDX → dataset derived by cutting at '_'
# ---------------------------------------------------------------------------
task SubsetChunk {
  input {
    String vcf
    String prefix
    String region
    File   resolved_mapping
    Int    chunk_mb
    Int    cpu     = 4
    Int    disk_gb = chunk_mb * 3 / 1024 + 5
  }

  String dataset = sub(prefix, "_.*", "")
  String out     = prefix + ".vcf.gz"

  command <<<
  set -euo pipefail

  # Filter resolved_mapping to this dataset
  awk -v ds="~{dataset}" '
    NR==1{next}
    $3==ds && $2!="NA" && $2!="AMBIGUOUS" && $2!="" {print $1"\t"$2}
  ' ~{resolved_mapping} > mapping.tsv

  export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

  cut -f1 mapping.tsv > sample_list.txt

  # Stream from GCS with fresh token
  bcftools view -r "~{region}" --samples-file sample_list.txt --force-samples --threads ~{cpu} -Oz -o subset.vcf.gz "~{vcf}"

  # Reheader locally — no GCS, no token risk
  bcftools reheader --samples mapping.tsv --output "~{out}" subset.vcf.gz
  rm subset.vcf.gz
  bcftools index -t --threads ~{cpu} "~{out}"

  echo "DONE ~{out} variants=$(bcftools index -n ~{out}) samples=$(bcftools query -l ~{out} | wc -l)"
  >>>

  output {
    File out_vcf = out
    File out_tbi = out + ".tbi"
  }

  runtime {
    cpu:   cpu
    disks: "local-disk ~{disk_gb} HDD"
  }
}


# ---------------------------------------------------------------------------
# Concat all chunks for one dataset+chrom into a single VCF.
# Receives all chunk paths as strings (no localization), greps to select
# matching chunks, sorts by chunk index, concatenates.
# ---------------------------------------------------------------------------
task ConcatChromVCF {
  input {
    Array[String] all_vcf_paths
    String        dataset
    String        vcf
    String        chrom
    String        suffix
    Int           chrom_size_mb
    Int           cpu     = 4
    Int           disk_gb = chrom_size_mb * 2 / 1024 + 10
  }

  String base = sub(basename(vcf, ".vcf.gz"), "\\.QC_ANNOTATED$", "")
  String out  = base + ".QC_ANNOTATED_" + suffix + "_" + chrom + ".vcf.gz"

  command <<<
  set -euo pipefail

  grep "~{dataset}_~{chrom}_" ~{write_lines(all_vcf_paths)} | sort > chunks.txt
  echo "~{dataset} ~{chrom}: $(wc -l < chunks.txt) chunks" >&2

  bcftools concat \
    --file-list chunks.txt \
    --output-type z \
    --threads ~{cpu} \
    --output ~{out}

  bcftools index -t --threads ~{cpu} ~{out}
  echo "DONE ~{out} variants=$(bcftools index -n ~{out})" >&2
  >>>

  output {
    File out_vcf = out
    File out_tbi = out + ".tbi"
  }

  runtime {
    cpu:   cpu
    disks: "local-disk ~{disk_gb} HDD"
  }
}
