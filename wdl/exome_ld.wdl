version 1.0

# Per-chromosome FG + exome plink conversion and merge.
#
# Phases:
#   1. BuildFgRegions  — single task; queries FG VCF stats (index only), splits
#                        per-chrom positions into chunk_mb-sized chunks, builds
#                        exome union sample list.
#   2. SubsetFgChunk   — scatter over all chunks; streams one FG region →
#                        sample-subsetted VCF.gz chunk.
#   3. ConcatFgChrom   — scatter over chroms; concatenates chunks → one VCF.gz + tbi.
#   4. VcfToPlink (FG) — scatter over chroms; subsetted VCF.gz → plink via FUSE.
#   5. VcfToPlink (ex) — flat scatter N_datasets × N_chroms; exome VCF → plink,
#                        --exclude FG bim for that chrom.
#   6. MergeChrom      — scatter over chroms; merges FG + all exome in one plink1 call.

workflow exome_ld {
  input {
    String               fg_vcf_template       # gs://…/finngen_R14_chrCHROM.vcf.gz
    File                 fg_pheno_file         # finngen_R14_minimum…txt.gz — for --update-sex
    Array[Array[String]] exome_vcf_pairs       # [["PREFIX", "gs://…/chrCHROM.vcf.gz"], …]
    Array[String]        chroms                # e.g. ["1","2",…,"22","23"]
    File                 positions_bim         # plink BIM with all chroms; col1=chrom col4=pos
    File                 credible_groups   # text file with one group_report gs:// path per line
    Int                  chunk_mb         = 1000
    String        plink_conv_args  = "--double-id --allow-extra-chr --split-par hg38 --vcf-half-call h"
    String        plink_merge_args = "--allow-extra-chr"
    String        out_prefix       = "finngen_R14_exome"
    String        ld_params        = "--ld-window-kb 1000 --ld-window-r2 0.05"
    File          annot                          # VEP annotation TSV.bgz
    String        filter_docker    = "eu.gcr.io/finngen-refinery-dev/exome_bioinf:ld.zstd"
    Int           mem_gb           = 36
    Int           cpu              = 8
    Int           disk_gb          = 200
  }

  # ── BuildLeads: generate leads table from group reports ─────────────────────
  Array[File] group_reports = read_lines(credible_groups)

  call BuildLeads {
    input:
      group_reports = group_reports,
      out_prefix    = out_prefix,
      docker        = filter_docker
  }

  # ── Phase 1: build chunk regions + exome union sample list ──────────────────
  call BuildFgRegions {
    input:
      fg_vcf_template      = fg_vcf_template,
      exome_vcf_pairs      = exome_vcf_pairs,
      positions_bim        = positions_bim,
      chroms               = chroms,
      chunk_mb             = chunk_mb
  }

  # ── Phase 2: scatter FG subset over chunks ───────────────────────────────────
  scatter (i in range(length(BuildFgRegions.tasks))) {
    call SubsetFgChunk {
      input:
        fg_vcf       = BuildFgRegions.tasks[i][0],
        chunk_prefix = BuildFgRegions.tasks[i][1],
        region       = BuildFgRegions.tasks[i][2],
        samples      = BuildFgRegions.union_samples,
        chunk_mb     = chunk_mb,
        cpu          = cpu
    }
  }

  # ── Phase 3+4: concat chunks then convert to plink, per chrom ───────────────
  scatter (ci in range(length(chroms))) {
    call ConcatFgChrom {
      input:
        all_chunk_vcfs = SubsetFgChunk.vcf,
        chrom          = chroms[ci],
        chrom_size_mb  = BuildFgRegions.chrom_sizes_mb[ci],
        cpu            = cpu
    }
    call VcfToPlink as FGtoPlink {
      input:
        vcf_template    = ConcatFgChrom.vcf,   # File coerced to String — no localisation
        chrom           = chroms[ci],
        out_prefix      = "fg_chr" + chroms[ci],
        fg_pheno_file   = fg_pheno_file,
        plink_conv_args = plink_conv_args,
        mem_gb          = mem_gb,
        cpu             = cpu,
        disk_gb         = ceil(size(ConcatFgChrom.vcf, "GB") * 2) + 5
    }
  }

  # ── Phase 5: exome VCFs → plink, exclude FG variants ────────────────────────
  scatter (i in range(length(exome_vcf_pairs) * length(chroms))) {
    Int vcf_i   = i / length(chroms)
    Int chrom_i = i % length(chroms)
    call VcfToPlink as ExomeToPlink {
      input:
        vcf_template    = exome_vcf_pairs[vcf_i][1],
        chrom           = chroms[chrom_i],
        out_prefix      = exome_vcf_pairs[vcf_i][0] + "_chr" + chroms[chrom_i],
        exclude_bim     = FGtoPlink.bim[chrom_i],
        fg_pheno_file   = fg_pheno_file,
        plink_conv_args = plink_conv_args,
        mem_gb          = mem_gb,
        cpu             = cpu,
        disk_gb         = disk_gb
    }
  }

  # ── Phase 6: merge FG + all exome per chrom ──────────────────────────────────
  scatter (ci in range(length(chroms))) {
    Float merge_fg_gb    = size(FGtoPlink.bed[ci], "GB")
    Float merge_exome_gb = size(ExomeToPlink.bed, "GB") / length(chroms)
    Float merge_total_gb = merge_fg_gb + merge_exome_gb

    call MergeChrom {
      input:
      chrom            = chroms[ci],
      out_prefix       = out_prefix,
      all_fg_beds      = FGtoPlink.bed,
      all_exome_beds   = ExomeToPlink.bed,
      plink_merge_args = plink_merge_args,
      mem_gb           = ceil(merge_total_gb * 2) + 4,
      cpu              = cpu,
      disk_gb          = ceil(merge_total_gb * 2) + 20
      }

      
    call PlinkToVcf {
      input:
      plink      = MergeChrom.plink,
      out_prefix = out_prefix,
      chrom      = chroms[ci],
      mem_gb     = mem_gb,
      cpu        = cpu,
      disk_gb    = ceil(merge_total_gb * 2) + 20
    }
      
    call ComputeLd {
      input:
      plink        = MergeChrom.plink,
      fg_bim       = FGtoPlink.bim[ci],   # File coerced to String — no localisation
      ld_params    = ld_params,
      chrom        = chroms[ci],
      cpu          = cpu*2,
      disk_gb      = if chroms[ci] == "6" then 100 else 50
    }

    call FilterLd {
      input:
      vcor_zst = ComputeLd.vcor,
      fg_bim   = FGtoPlink.bim[ci],
      annot    = annot,
      afreq    = MergeChrom.afreq,
      chrom    = chroms[ci],
      docker   = filter_docker
    }

    call AnnotateLd {
      input:
      leads_tsv = BuildLeads.leads_tsv,
      ld_file   = FilterLd.ld,
      docker    = filter_docker
    }

  }

  call GatherAnnotatedLd {
    input:
      annotated_chroms = AnnotateLd.annotated,
      ld_chroms        = FilterLd.ld,
      out_prefix       = out_prefix,
      docker           = filter_docker
  }

  output {
    Array[Array[File]] merged_plink    = MergeChrom.plink
    Array[File]        afreq           = MergeChrom.afreq
    Array[File]        ld_results      = FilterLd.ld
    Array[File]        merged_vcf      = PlinkToVcf.vcf
    File               leads_tsv_out   = BuildLeads.leads_tsv
    File               annotated_ld    = GatherAnnotatedLd.annotated
    File               raw_ld          = GatherAnnotatedLd.raw_ld
  }
}


# ---------------------------------------------------------------------------
# Single task: query FG VCF stats (index only), split per-chrom positions
# into size-based chunks, build exome union sample list.
#
# Output tasks.tsv columns: fg_vcf_gs_path | chunk_prefix | region
# region format: chr{chrom}:start-end  (FinnGen VCFs use chr-prefixed names)
# chrom_sizes_mb: expected output MB per chrom after sample_reduction.
# ---------------------------------------------------------------------------
task BuildFgRegions {
  input {
    String               fg_vcf_template
    Array[Array[String]] exome_vcf_pairs
    File                 positions_bim
    Array[String]        chroms
    Int                  chunk_mb
  }

  command <<<
  set -euo pipefail

  # ── inputs ──────────────────────────────────────────────────────────────────
  FG_VCF_TEMPLATE="~{fg_vcf_template}"
  EXOME_PAIRS_TSV="~{write_tsv(exome_vcf_pairs)}"
  CHROMS_FILE="~{write_lines(chroms)}"
  POSITIONS_BIM="~{positions_bim}"
  CHUNK_MB="~{chunk_mb}"

  export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
  export HTS_HTTP_VERSION=1.1

  # ── exome union sample list (first chrom to resolve template) ───────────────
  first_chrom=$(head -1 "$CHROMS_FILE")
  while IFS=$'\t' read -r dataset vcf_tmpl; do
    bcftools query -l "${vcf_tmpl/CHROM/$first_chrom}"
  done < "$EXOME_PAIRS_TSV" | sort -u > union_samples.txt
  n_exome=$(wc -l < union_samples.txt)
  echo "Exome union: $n_exome samples" >&2

  # ── FG sample count + bpv from first chrom VCF ──────────────────────────────
  first_vcf="${FG_VCF_TEMPLATE/CHROM/$first_chrom}"
  n_fg=$(bcftools query -l "$first_vcf" | wc -l)
  echo "FG samples: $n_fg" >&2

  sample_reduction=$(awk -v ne="$n_exome" -v nf="$n_fg" 'BEGIN{printf "%.6f", ne/nf}')
  echo "sample_reduction=$sample_reduction" >&2

  file_bytes=$(gsutil du "$first_vcf" | awk '{print $1}')
  n_total=$(bcftools index --stats "$first_vcf" | awk '{sum+=$3} END{print sum}')
  bpv=$(awk -v fb="$file_bytes" -v nt="$n_total" 'BEGIN{printf "%.4f", fb/nt}')
  echo "bpv=$bpv ($(( file_bytes/1024/1024 ))MB / $n_total variants from chr$first_chrom)" >&2

  n_per_chunk=$(awk -v bpv="$bpv" -v sr="$sample_reduction" -v mb="$CHUNK_MB" \
    'BEGIN{n=int(mb*1024*1024/(bpv*sr)); print (n<1?1:n)}')
  echo "n_per_chunk=$n_per_chunk" >&2

  # ── split BIM by chrom in one pass; derive VCF chrom name from variant ID ───
  awk '{split($2,a,"_"); print $4 > "pos_"$1".txt"; vcf_chrom[$1]=a[1]}
       END{for(c in vcf_chrom) print c"\t"vcf_chrom[c] > "chrom_map.txt"}' "$POSITIONS_BIM"

  # ── per-chrom: size estimate + chunk regions ─────────────────────────────────
  while IFS= read -r chrom; do
    vcf="${FG_VCF_TEMPLATE/CHROM/$chrom}"
    pos_file="pos_${chrom}.txt"
    vcf_chrom=$(awk -v c="$chrom" '$1==c{print $2; exit}' chrom_map.txt)

    n_pos=$(wc -l < "$pos_file")
    chrom_mb=$(awk -v bpv="$bpv" -v n="$n_pos" -v sr="$sample_reduction" \
      'BEGIN{printf "%.0f", bpv*n*sr/1024/1024}')
    echo "$chrom_mb" >> chrom_sizes_mb.txt

    contig_end=$(bcftools view -h "$vcf" | grep "##contig=<ID=${vcf_chrom}," | grep -o 'length=[0-9]*' | cut -d= -f2 || true)

    split -l "$n_per_chunk" "$pos_file" "split_${chrom}_"

    mapfile -t split_files < <(ls "split_${chrom}_"* | sort)
    n_files=${#split_files[@]}
    chunk=0; start=1
    for i in "${!split_files[@]}"; do
      f="${split_files[$i]}"
      end=$(tail -n1 "$f")
      [[ $i -eq $((n_files - 1)) && -n "$contig_end" ]] && end=$contig_end
      prefix="fg_${chrom}_$(printf '%03d' $chunk)"
      printf "%s\t%s\t%s:%d-%d\n" "$vcf" "$prefix" "$vcf_chrom" "$start" "$end"
      start=$((end + 1))
      chunk=$((chunk + 1))
    done
    rm "$pos_file" "split_${chrom}_"*
  done < "$CHROMS_FILE" > tasks.tsv

  echo "Total tasks: $(wc -l < tasks.tsv)" >&2
  >>>

  output {
    Array[Array[String]] tasks          = read_tsv("tasks.tsv")
    Array[Int]           chrom_sizes_mb = read_lines("chrom_sizes_mb.txt")
    File                 union_samples  = "union_samples.txt"
  }

}


# ---------------------------------------------------------------------------
# Stream one FG VCF region directly from GCS, subset to exome union samples.
# Retries up to 3× on transient GCS errors.
# ---------------------------------------------------------------------------
task SubsetFgChunk {
  input {
    String fg_vcf
    String chunk_prefix
    String region
    File   samples
    Int    chunk_mb
    Int    cpu     = 4
    Int    disk_gb = chunk_mb * 2 / 1024 + 5
  }

  String out = chunk_prefix + ".vcf.gz"

  command <<<
  set -euo pipefail

  fuse_vcf=$(echo "~{fg_vcf}" | sed 's|gs://[^/]*/|/mnt/disks/gcs/|')

  bcftools view \
    -r "~{region}" -t "~{region}" \
    -S ~{samples} --force-samples \
    --threads ~{cpu} \
    -Oz -o "~{out}" \
    "$fuse_vcf"
  echo "DONE ~{out}: $(bcftools index -n ~{out} 2>/dev/null || echo '?') variants" >&2
  >>>

  output {
    File vcf = out
  }

  runtime {
    cpu:   cpu
    disks: "local-disk ~{disk_gb} HDD"
  }
}


# ---------------------------------------------------------------------------
# Concat all FG chunks for one chrom into a single indexed VCF.gz.
# all_chunk_vcfs is Array[File] coerced to String — no localisation.
# ---------------------------------------------------------------------------
task ConcatFgChrom {
  input {
    Array[String] all_chunk_vcfs
    String        chrom
    Int           chrom_size_mb
    Int           cpu     = 8
    Int           disk_gb = chrom_size_mb * 2 / 1024 + 10
  }

  String out = "fg_" + chrom + "_subset.vcf.gz"

  command <<<
  set -euo pipefail
  export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

  grep "/fg_~{chrom}_" ~{write_lines(all_chunk_vcfs)} | \
    awk -F/ '{print $NF "\t" $0}' | sort -k1,1 | cut -f2- > chunks.txt
  echo "~{chrom}: $(wc -l < chunks.txt) chunks" >&2

  bcftools concat -n \
    --file-list chunks.txt \
    --output-type z \
    --output "~{out}"

  bcftools index -t --threads ~{cpu} "~{out}"
  echo "DONE ~{out}: $(bcftools index -n ~{out}) variants" >&2
  >>>

  output {
    File vcf = out
    File tbi = out + ".tbi"
  }

  runtime {
    cpu:   cpu
    disks: "local-disk ~{disk_gb} HDD"
  }
}


# ---------------------------------------------------------------------------
# VCF/VCF.gz → plink via GCS FUSE.
# vcf_template may carry a CHROM placeholder (exome) or be an already-resolved
# gs:// path (FG subsetted VCF, File coerced to String at workflow level).
# ---------------------------------------------------------------------------
task VcfToPlink {
  input {
    String        vcf_template     # gs://… with CHROM placeholder, or resolved gs:// path
    String        chrom
    String        out_prefix
    File?         exclude_bim
    File          fg_pheno_file
    String        plink_conv_args
    Int           mem_gb
    Int           cpu
    Int           disk_gb
  }

  Int plink_mem_mb = if mem_gb > 6 then (mem_gb - 4) * 1024 else 2048

  command <<<
  set -euo pipefail
  export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

  resolve_fuse() {
    local fuse
    fuse=$(echo "$1" | sed "s/CHROM/~{chrom}/; s|gs://[^/]*/|/mnt/disks/gcs/|")
    [[ -f "$fuse" ]] || { echo "ERROR: not found: $1 (chrom ~{chrom})" >&2; exit 1; }
    echo "$fuse"
  }

  fuse_path=$(resolve_fuse "~{vcf_template}")
  echo "FUSE path: $fuse_path" >&2

  zcat -f ~{fg_pheno_file} | awk -F'\t' 'NR>1{sex=($5=="male")?1:($5=="female")?2:0; print $1,$1,sex}' > sex_update.txt

  plink2 \
    --vcf        "$fuse_path" \
    ~{if defined(exclude_bim) then "--exclude " + select_first([exclude_bim]) else ""} \
    --update-sex sex_update.txt \
    --make-bed \
    --memory     ~{plink_mem_mb} \
    --threads    ~{cpu} \
    --out        "~{out_prefix}" \
    ~{plink_conv_args}

  sed -i 's/\tchrX_/\tchr23_/g' "~{out_prefix}.bim"
  echo "~{out_prefix}: $(wc -l < ~{out_prefix}.bim) variants, $(wc -l < ~{out_prefix}.fam) samples" >&2
  >>>

  output {
    File bed = out_prefix + ".bed"
    File bim = out_prefix + ".bim"
    File fam = out_prefix + ".fam"
    File log = out_prefix + ".log"
  }

  runtime {
    cpu:    cpu
    memory: mem_gb + " GB"
    disks:  "local-disk ~{disk_gb} HDD"
  }
}


# ---------------------------------------------------------------------------
# Merge FG + all exome plink files for one chrom via plink1.
# Inputs are String (no localisation); files are reached via GCS FUSE.
# ---------------------------------------------------------------------------
task MergeChrom {
  input {
    String        chrom
    String        out_prefix
    Array[String] all_fg_beds     # Array[File] coerced to String — no localisation
    Array[String] all_exome_beds
    String        plink_merge_args
    Int           mem_gb
    Int           cpu
    Int           disk_gb
  }

  Int    plink_mem_mb = if mem_gb > 6 then (mem_gb - 4) * 1024 else 2048
  String out          = out_prefix + "_chr" + chrom

  command <<<
  set -euo pipefail

  to_fuse() { echo "$1" | sed 's|gs://[^/]*/|/mnt/disks/gcs/|'; }

  fg_gcs=$(grep "chr~{chrom}\." ~{write_lines(all_fg_beds)})
  fg_fuse=$(to_fuse "$fg_gcs")
  fg_prefix="${fg_fuse%.bed}"
  echo "FG: $fg_prefix" >&2

  grep "chr~{chrom}\." ~{write_lines(all_exome_beds)} | while read gcs; do
    fuse=$(to_fuse "$gcs")
    prefix="${fuse%.bed}"
    echo "${prefix}.bed ${prefix}.bim ${prefix}.fam"
  done > merge_list.txt
  echo "Exome filesets: $(wc -l < merge_list.txt)" >&2

  plink \
    --bfile      "$fg_prefix" \
    --merge-list merge_list.txt \
    --make-bed \
    --memory     ~{plink_mem_mb} \
    --threads    ~{cpu} \
    --out        "~{out}" \
    ~{plink_merge_args}

  echo "~{out}: $(wc -l < ~{out}.bim) variants, $(wc -l < ~{out}.fam) samples" >&2

  plink2 \
    --bfile   "~{out}" \
    --freq \
    --memory  ~{plink_mem_mb} \
    --threads ~{cpu} \
    --out     "~{out}"
  >>>

  output {
    Array[File] plink = [out + ".bed", out + ".bim", out + ".fam", out + ".log"]
    File        afreq = out + ".afreq"
  }

  runtime {
    cpu:    cpu
    memory: mem_gb + " GB"
    disks:  "local-disk ~{disk_gb} HDD"
  }
}


# ---------------------------------------------------------------------------
# Compute LD between exome and FG variants in the merged plink fileset.
# Annotates the merged BIM with fg/exome prefixes so plink2 output can be
# filtered to cross-dataset pairs only.
# ---------------------------------------------------------------------------
task ComputeLd {
  input {
    Array[String] plink      # [bed, bim, fam, log] from MergeChrom — no localisation
    String        fg_bim
    String        ld_params
    String        chrom
    Int           cpu     = 8
    Int           mem_gb  = 16
    Int           disk_gb = 20
  }

  Int    plink_mem_mb = if mem_gb > 6 then (mem_gb - 4) * 1024 else 2048
  String out_root     = "exome_finngen_ld_" + chrom

  command <<<
  set -euo pipefail

  to_fuse() { echo "$1" | sed 's|gs://[^/]*/|/mnt/disks/gcs/|'; }

  bed=$(to_fuse "~{plink[0]}")
  bim=$(to_fuse "~{plink[1]}")
  fam=$(to_fuse "~{plink[2]}")
  fg_bim=$(to_fuse "~{fg_bim}")

  plink2 \
    --bed "$bed" \
    --bim "$bim" \
    --fam "$fam" \
    --ld-snp-list "$fg_bim" \
    --r2-unphased zs ~{ld_params} \
    --memory ~{plink_mem_mb} \
    --threads ~{cpu} \
    --zst-level 1 \
    --out ld

  mv ld.vcor.zst ~{out_root}.vcor.zst
  >>>


  output {
    File vcor = out_root + ".vcor.zst"
  }

  runtime {
    cpu:    cpu
    memory: mem_gb + " GB"
    disks:  "local-disk ~{disk_gb} HDD"
  }
}


# ---------------------------------------------------------------------------
# Filter raw plink2 .vcor to FG→exome pairs and flag coding status.
# ---------------------------------------------------------------------------
task FilterLd {
  input {
    File   vcor_zst
    File   fg_bim
    File   annot
    File   afreq
    String chrom
    String docker
    Int    mem_gb  = 16
    Int    disk_gb = 20
  }

  String out = "exome_finngen_ld_" + chrom + ".ld.tsv.gz"

  command <<<
  set -euo pipefail

  zstd -d ~{vcor_zst} -c | bgzip > tmp.vcor.gz
  echo "~{chrom}: $(zcat tmp.vcor.gz | tail -n+2 | wc -l) raw pairs" >&2

  python3 << 'PYEOF'
  import sys, pandas as pd

  CODING = {"missense_variant","stop_gained","frameshift_variant","splice_acceptor_variant",
            "splice_donor_variant","start_lost","stop_lost","inframe_insertion","inframe_deletion"}

  fg_ids       = {line.split()[1] for line in open("~{fg_bim}")}
  annot        = pd.read_csv("~{annot}", sep='\t', usecols=['rsid','most_severe','gene_most_severe'], compression='gzip')
  annot['rsid'] = annot['rsid'].str.replace('chrX_', 'chr23_', regex=False)
  conseq_map   = annot[annot['most_severe'].isin(CODING)].drop_duplicates('rsid').set_index('rsid')['most_severe'].to_dict()
  gene_map     = annot.drop_duplicates('rsid').set_index('rsid')['gene_most_severe'].to_dict()
  afreq        = pd.read_csv("~{afreq}", sep='\t')
  afreq.columns = afreq.columns.str.lstrip('#')
  af_map       = afreq.set_index('ID')['ALT_FREQS'].to_dict()
  print(f"~{chrom}: {len(fg_ids)} FG IDs, {len(conseq_map)} coding, {len(af_map)} AF entries", file=sys.stderr)

  import gzip
  with gzip.open("~{out}", 'wt') as out:
    for i, chunk in enumerate(pd.read_csv("tmp.vcor.gz", sep='\t', chunksize=500_000, usecols=['ID_A','ID_B','UNPHASED_R2'])):
      df = chunk[~chunk['ID_B'].isin(fg_ids)].rename(columns={'ID_A':'FG_SNP','ID_B':'EXOME_SNP','UNPHASED_R2':'R2'})
      df['exome_consequence']  = df['EXOME_SNP'].map(conseq_map).fillna('NA')
      df['EXOME_AF']           = df['EXOME_SNP'].map(af_map)
      df['exome_nearest_gene'] = df['EXOME_SNP'].map(gene_map)
      df.to_csv(out, sep='\t', index=False, header=(i==0))
  print(f"~{chrom}: done", file=sys.stderr)
  PYEOF

  rm tmp.vcor.gz
  >>>

  output {
    File ld = out
  }

  runtime {
    docker: docker
    memory: mem_gb + " GB"
    disks:  "local-disk ~{disk_gb} HDD"
  }
}


# ---------------------------------------------------------------------------
# Convert merged plink fileset to bgzipped VCF.
# plink is Array[String] (no localisation); files are reached via GCS FUSE.
# ---------------------------------------------------------------------------
task PlinkToVcf {
  input {
    Array[String] plink       # [bed, bim, fam, log] from MergeChrom — no localisation
    String        out_prefix
    String        chrom
    Int           mem_gb  = 16
    Int           cpu     = 4
    Int           disk_gb = 50
  }

  Int    plink_mem_mb = if mem_gb > 6 then (mem_gb - 4) * 1024 else 2048
  String out          = out_prefix + "_chr" + chrom

  command <<<
  set -euo pipefail

  to_fuse() { echo "$1" | sed 's|gs://[^/]*/|/mnt/disks/gcs/|'; }
  bed=$(to_fuse "~{plink[0]}")
  bfile="${bed%.bed}"

  plink2 \
    --bfile      "$bfile" \
    --export     vcf bgz \
    --output-chr chrM \
    --memory     ~{plink_mem_mb} \
    --threads    ~{cpu} \
    --out        "~{out}"

  echo "~{out}: done" >&2
  >>>

  output {
    File vcf = out + ".vcf.gz"
  }

  runtime {
    cpu:    cpu
    memory: mem_gb + " GB"
    disks:  "local-disk ~{disk_gb} HDD"
  }
}


# ---------------------------------------------------------------------------
# Build leads table from FinnGen group reports.
# ---------------------------------------------------------------------------
task BuildLeads {
  input {
    Array[File] group_reports
    String      out_prefix
    String      docker
    Int         mem_gb  = 16
    Int         disk_gb = 50
  }

  command <<<
  set -euo pipefail

  python3 << 'PYEOF'
  import sys, pandas as pd

  READ_COLS = ["phenotype_abbreviation", "locus_id", "lead_mlogp", "lead_beta",
               "lead_af_alt", "good_cs", "best_coding_var", "functional_variants_relaxed"]
  KEEP_COLS = ["PHENO", "locus_id", "lead_mlogp", "lead_beta", "lead_af_alt", "good_cs"]
  paths = open("~{write_lines(group_reports)}").read().splitlines()
  print(f"{len(paths)} group_report files", file=sys.stderr)

  dfs = []
  for path in paths:
    try:
      dfs.append(pd.read_csv(path, sep="\t", usecols=READ_COLS, low_memory=False))
    except Exception as e:
      print(f"WARNING: could not read {path}: {e}", file=sys.stderr)
  df = pd.concat(dfs, ignore_index=True)
  df = df.rename(columns={"phenotype_abbreviation": "PHENO"})

  has_coding     = df["best_coding_var"].notna()
  has_functional = df["functional_variants_relaxed"].notna()
  df["cs_type"] = "NA"
  df.loc[has_functional & ~has_coding, "cs_type"] = "functional_relaxed"
  df.loc[has_coding,                   "cs_type"] = "coding"

  def parse(val):
    if pd.isna(val): return None, None
    p = val.split(";")[0].split("|")
    return p[0], float(p[3])
  parsed = pd.DataFrame(df["functional_variants_relaxed"].map(parse).tolist(),
                        index=df.index, columns=["functional_var", "functional_var_r2"])
  df["functional_var"]    = parsed["functional_var"]
  df["functional_var_r2"] = parsed["functional_var_r2"]

  leads = df[KEEP_COLS + ["cs_type", "functional_var", "functional_var_r2"]]
  leads.to_csv("~{out_prefix}.leads.tsv.gz", sep="\t", index=False, compression="gzip", na_rep="NA")
  print(f"{len(leads)} CSs across {leads['PHENO'].nunique()} phenotypes", file=sys.stderr)
  PYEOF
  >>>

  output {
    File leads_tsv = out_prefix + ".leads.tsv.gz"
  }

  runtime {
    docker: docker
    memory: mem_gb + " GB"
    disks:  "local-disk ~{disk_gb} HDD"
  }
}


# ---------------------------------------------------------------------------
# Annotate per-chrom LD file by joining to the leads table.
# ---------------------------------------------------------------------------
task AnnotateLd {
  input {
    File   leads_tsv
    File   ld_file
    String docker
    Int    mem_gb  = 16
    Int    disk_gb = 50
  }

  String out = sub(basename(ld_file), "\\.ld\\.tsv\\.gz$", ".ld_annotated.tsv.gz")

  command <<<
  set -euo pipefail

  python3 << 'PYEOF'
  import gzip, sys, pandas as pd
  from tqdm import tqdm

  leads = pd.read_csv("~{leads_tsv}", sep="\t")
  print(f"{len(leads)} CSs loaded", file=sys.stderr)

  n_rows = 0
  header_written = False
  with gzip.open("~{out}", "wt") as out:
    for chunk in tqdm(pd.read_csv("~{ld_file}", sep="\t", chunksize=500_000),
                      desc="annotating", unit="chunk", file=sys.stderr):
      merged = chunk.merge(leads, left_on="FG_SNP", right_on="locus_id", how="inner").drop(columns="locus_id")
      if not merged.empty:
        merged.to_csv(out, sep="\t", index=False, header=not header_written, na_rep="NA")
        header_written = True
        n_rows += len(merged)
  print(f"{n_rows} annotated rows written to ~{out}", file=sys.stderr)
  PYEOF
  >>>

  output {
    File annotated = out
  }

  runtime {
    docker: docker
    memory: mem_gb + " GB"
    disks:  "local-disk ~{disk_gb} HDD"
  }
}


# ---------------------------------------------------------------------------
# Gather per-chrom annotated LD files into one.
# ---------------------------------------------------------------------------
task GatherAnnotatedLd {
  input {
    Array[String] annotated_chroms
    Array[String] ld_chroms
    String        out_prefix
    String        docker
    Int           mem_gb  = 8
  }

  String out     = out_prefix + ".ld_annotated.tsv.gz"
  String out_raw = out_prefix + ".ld.tsv.gz"

  command <<<
  set -euo pipefail

  sed 's|gs://[^/]*/|/mnt/disks/gcs/|' ~{write_lines(annotated_chroms)} > annotated_fuse.txt
  sed 's|gs://[^/]*/|/mnt/disks/gcs/|' ~{write_lines(ld_chroms)}        > ld_fuse.txt

  mapfile -t annotated_files < annotated_fuse.txt
  mapfile -t ld_files < ld_fuse.txt

  echo "gathering ${#annotated_files[@]} annotated files" >&2
  zcat "${annotated_files[0]}" | bgzip -c > ~{out}
  for f in "${annotated_files[@]:1}"; do
    echo "annotated: $f" >&2
    zcat "$f" | sed -E 1d | bgzip -c >> ~{out}
  done
  echo "wrote ~{out}" >&2

  echo "gathering ${#ld_files[@]} raw ld files" >&2
  zcat "${ld_files[0]}" | bgzip -c > ~{out_raw}
  for f in "${ld_files[@]:1}"; do
    echo "raw ld: $f" >&2
    zcat "$f" | sed -E 1d | bgzip -c >> ~{out_raw}
  done
  echo "wrote ~{out_raw}" >&2
  >>>

  output {
    File annotated = out
    File raw_ld    = out_raw
  }

  runtime {
    docker: docker
    memory: mem_gb + " GB"
    disks:  "local-disk 100 HDD"
  }
}
