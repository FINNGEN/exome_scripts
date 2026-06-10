version 1.0

# Per-chromosome FG + exome plink conversion and merge via GCS FUSE.
# Mount: gs://BUCKET/path → /mnt/disks/gcs/path
#
# Phases:
#   1. SubsetFgVcf  — scatter over chroms; FG VCF → sample-subsetted VCF.gz via FUSE.
#                     Resolves chr23/chrX by probing FUSE paths (templates passed as-is).
#   2. VcfToPlink   — scatter over chroms; subsetted VCF.gz → plink. Disk sized from file.
#   3. VcfToPlink   — flat scatter over N_datasets × N_chroms; exome VCF → plink,
#                     --exclude FG bim for that chrom. Same task, same resolve logic.
#   4. MergeChrom   — scatter over chroms; merges FG + all exome in one plink1 call.

workflow exome_ld {
  input {
    String               fg_vcf_template   # gs://…/finngen_R14_chrCHROM.vcf.gz
    File                 fg_pheno_file     # finngen_R14_minimum…txt.gz — for --update-sex
    Array[Array[String]] exome_vcf_pairs  # [["PREFIX", "gs://…/chrCHROM.vcf.gz"], …]
    Array[String]        chroms           # e.g. ["1","2",…,"22","23"]
    String        plink_conv_args  = "--double-id --allow-extra-chr --split-par hg38 --vcf-half-call h"
    String        plink_merge_args = "--allow-extra-chr"
    Int           plink_mem_gb     = 32
    Int           cpu              = 8
    Int           disk_gb          = 200
  }

  # ── Phase 1+2: subset FG VCF to exome samples, convert to plink ─────────────
  scatter (chrom in chroms) {
    call SubsetFgVcf {
      input:
        fg_vcf_template     = fg_vcf_template,
        exome_vcf_templates = transpose(exome_vcf_pairs)[1],
        chrom               = chrom,
        cpu                 = cpu,
        disk_gb             = disk_gb
    }
    call VcfToPlink as FGtoPlink {
      input:
        vcf_template    = SubsetFgVcf.vcf,   # File coerced to String — no localisation
        chrom           = chrom,
        out_prefix      = "fg_chr" + chrom,
        fg_pheno_file   = fg_pheno_file,
        plink_conv_args = plink_conv_args,
        plink_mem_gb    = plink_mem_gb,
        cpu             = cpu,
        disk_gb         = ceil(size(SubsetFgVcf.vcf, "GB") * 2) + 5
    }
  }

  # ── Phase 3: exome VCFs → plink, exclude FG variants ───────────────────────
  scatter (i in range(length(exome_vcf_pairs) * length(chroms))) {
    Int vcf_i   = i / length(chroms)
    Int chrom_i = i % length(chroms)
    call VcfToPlink as ExomeToPlink {
      input:
        vcf_template    = exome_vcf_pairs[vcf_i][1],   # CHROM resolved inside task
        chrom           = chroms[chrom_i],
        out_prefix      = exome_vcf_pairs[vcf_i][0] + "_chr" + chroms[chrom_i],
        exclude_bim     = FGtoPlink.bim[chrom_i],
        fg_pheno_file   = fg_pheno_file,
        plink_conv_args = plink_conv_args,
        plink_mem_gb    = plink_mem_gb,
        cpu             = cpu,
        disk_gb         = disk_gb
    }
  }

  # ── Phase 4: merge FG + all exome per chrom ─────────────────────────────────
  scatter (chrom in chroms) {
    call MergeChrom {
      input:
        chrom            = chrom,
        all_fg_beds      = FGtoPlink.bed,
        all_exome_beds   = ExomeToPlink.bed,
        plink_merge_args = plink_merge_args,
        plink_mem_gb     = plink_mem_gb,
        cpu              = cpu,
        disk_gb          = disk_gb
    }
  }

  output {
    Array[File] merged_beds = MergeChrom.bed
    Array[File] merged_bims = MergeChrom.bim
    Array[File] merged_fams = MergeChrom.fam
  }
}


# ---------------------------------------------------------------------------
# Subset FG VCF to exome-union samples via GCS FUSE.
# Both fg_vcf_template and exome_vcf_templates carry a CHROM placeholder;
# resolve_fuse substitutes the chrom and falls back to 23↔X if the file
# isn't found at the first substitution.
# ---------------------------------------------------------------------------
task SubsetFgVcf {
  input {
    String        fg_vcf_template
    Array[String] exome_vcf_templates
    String        chrom
    Int           cpu
    Int           disk_gb
  }

  command <<<
  set -euo pipefail
  export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

  # Substitute CHROM, map to FUSE, fall back to 23↔X if not found.
  resolve_fuse() {
    local fuse
    fuse=$(echo "$1" | sed "s/CHROM/~{chrom}/; s|gs://[^/]*/|/mnt/disks/gcs/|")
    [[ -f "$fuse" ]] && { echo "$fuse"; return 0; }
    if [[ "~{chrom}" == "23" ]]; then
      fuse=$(echo "$1" | sed "s/CHROM/X/; s|gs://[^/]*/|/mnt/disks/gcs/|")
    elif [[ "~{chrom}" == "X" ]]; then
      fuse=$(echo "$1" | sed "s/CHROM/23/; s|gs://[^/]*/|/mnt/disks/gcs/|")
    fi
    [[ -f "$fuse" ]] || { echo "WARN: not found: $1 (chrom ~{chrom})" >&2; return 1; }
    echo "$fuse"
  }

  # Build exome union sample list
  > vcf_list.txt
  while IFS= read -r tmpl; do
    fuse=$(resolve_fuse "$tmpl") && echo "$fuse" >> vcf_list.txt || true
  done < ~{write_lines(exome_vcf_templates)}
  echo "found $(wc -l < vcf_list.txt) exome VCFs for chr~{chrom}" >&2

  while IFS= read -r vcf; do bcftools query -l "$vcf"; done < vcf_list.txt | sort -u > union_samples.txt
  echo "exome union: $(wc -l < union_samples.txt) samples" >&2

  # Subset FG VCF to exome samples
  fg_fuse=$(resolve_fuse "~{fg_vcf_template}")
  echo "subsetting FG: $fg_fuse" >&2
  bcftools view -S union_samples.txt --force-samples -O z -o "fg_chr~{chrom}_subset.vcf.gz" "$fg_fuse"
  >>>

  output {
    File vcf = "fg_chr~{chrom}_subset.vcf.gz"
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
# resolve_fuse handles both cases plus the 23↔X fallback.
# ---------------------------------------------------------------------------
task VcfToPlink {
  input {
    String        vcf_template     # gs://… with CHROM placeholder, or resolved gs:// path
    String        chrom
    String        out_prefix
    File?         exclude_bim
    File          fg_pheno_file
    String        plink_conv_args
    Int           plink_mem_gb
    Int           cpu
    Int           disk_gb
  }

  Int plink_mem_mb = plink_mem_gb * 1024

  command <<<
  set -euo pipefail
  export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

  resolve_fuse() {
    local fuse
    fuse=$(echo "$1" | sed "s/CHROM/~{chrom}/; s|gs://[^/]*/|/mnt/disks/gcs/|")
    [[ -f "$fuse" ]] && { echo "$fuse"; return 0; }
    if [[ "~{chrom}" == "23" ]]; then
      fuse=$(echo "$1" | sed "s/CHROM/X/; s|gs://[^/]*/|/mnt/disks/gcs/|")
    elif [[ "~{chrom}" == "X" ]]; then
      fuse=$(echo "$1" | sed "s/CHROM/23/; s|gs://[^/]*/|/mnt/disks/gcs/|")
    fi
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

  echo "~{out_prefix}: $(wc -l < ~{out_prefix}.bim) variants, $(wc -l < ~{out_prefix}.fam) samples" >&2
  >>>

  output {
    File bed = out_prefix + ".bed"
    File bim = out_prefix + ".bim"
    File fam = out_prefix + ".fam"
  }

  runtime {
    cpu:   cpu
    disks: "local-disk ~{disk_gb} HDD"
  }
}


# ---------------------------------------------------------------------------
# Merge FG + all exome plink files for one chrom via plink1.
# Inputs are String (no localisation); files are reached via GCS FUSE.
# ---------------------------------------------------------------------------
task MergeChrom {
  input {
    String        chrom
    Array[String] all_fg_beds     # Array[File] coerced to String — no localisation
    Array[String] all_exome_beds
    String        plink_merge_args
    Int           plink_mem_gb
    Int           cpu
    Int           disk_gb
  }

  Int plink_mem_mb = plink_mem_gb * 1024

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
    --out        "merged_chr~{chrom}" \
    ~{plink_merge_args}

  echo "merged_chr~{chrom}: $(wc -l < merged_chr~{chrom}.bim) variants, $(wc -l < merged_chr~{chrom}.fam) samples" >&2
  >>>

  output {
    File bed = "merged_chr~{chrom}.bed"
    File bim = "merged_chr~{chrom}.bim"
    File fam = "merged_chr~{chrom}.fam"
  }

  runtime {
    cpu:   cpu
    disks: "local-disk ~{disk_gb} HDD"
  }
}
