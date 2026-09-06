version 1.0

# Single self-contained workflow: kinship-based ID resolution (mandatory) plus
# an optional subset+rename stage — formerly two separate WDLs
# (exome_duplicates.wdl, exome_rename.wdl), now fully inlined into one file.
#
# Run with run_rename = false first to compute and inspect the ID mapping
# (id_mapping / id_mapping_stats / id_mapping_flowchart) without paying the
# cost of the rename stage. Once the mapping looks right, resubmit the same
# workflow with only run_rename = true changed — Cromwell call-caches every
# upstream call (none of their inputs changed), so only the newly enabled
# rename stage actually runs.
#
# Mandatory stage shape (formerly exome_duplicates.wdl):
#   MakeRegionSnplists
#         │
#         └─ FilterVCF ×(n_vcfs × n_regions)   [parallel scatter]
#                │
#          ConcatVCF                            [per VCF]
#                │
#          SubsetQuery (VCF → plink, shared HM3 variants)
#                │
#           ┌────┴────────────────┐
#      SubsetRef              FilterSNPs
#      (ref plink →           (QC filter query plink,
#       shared variants)       build final SNP list)
#           │                      │
#           └────────┬─────────────┘
#               PrepQuery / PrepRef
#               (subset to QC snplist, het-filter,
#                tag IDs, split into chunks)
#                      │
#                 KingShards
#                 (KING --duplicate, all chunk pairs)
#                      │
#               SummarizeKing
#               (per-sample duplicate summary)
#                      │
#                GatherResults
#                (combine summaries, plots, global stats,
#                 resolve ambiguities, final mapping + stats)
#
# Optional rename stage shape (formerly exome_rename.wdl), gated by run_rename:
#   QueryChromPositions[C]  one task per chrom — queries positions for all
#                           datasets in parallel (cached per chrom)
#   BuildAllRegions         single task — uses position files + file stats to
#                           build per-dataset size-based chunk TSV
#   SubsetChunk[K]          one task per (dataset × chrom × chunk)
#   ConcatChromVCF[D×C]     cross-product scatter — greps flat VCF path array
#                           to concat chunks for one dataset+chrom

workflow exome_reassign_ids {
  input {
    # ---- shared: same [DATASET, QC_ANNOTATED.vcf.gz] pairs feed both stages ----
    Array[Array[String]] vcf_pairs

    # ---- kinship / ID resolution (mandatory) ----
    File   plink_bed
    String plink_prefix
    Int    n_regions     = 100
    File?  berisa_blocks
    File?  aliases
    Int    target_snps   = 10000
    Float  max_het_F     = 0.3
    Int    chunk_size    = 10000
    String out_prefix    = "finngen_R14_exome"

    # ---- rename stage: subset + rename to FinnGen IDs (optional) ----
    Boolean       run_rename = false   # flip to true and resubmit once id_mapping looks right
    Array[String] chroms = [
      "chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
      "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
      "chr20","chr21","chr22","chrX","chrY"
    ]
    Int    chunk_mb = 500
    String suffix   = "fg_ids"
  }

  File plink_bim = sub(plink_bed, "\\.bed$", ".bim")
  File plink_fam = sub(plink_bed, "\\.bed$", ".fam")

  # =====================================================================
  # MANDATORY: kinship / duplicate detection
  # =====================================================================

  # Step 1: bin HM3 SNPs into N LD-block-merged regions
  call MakeRegionSnplists {
    input:
      bim           = plink_bim,
      n_regions     = n_regions,
      berisa_blocks = berisa_blocks
  }

  # Step 2: extract each (VCF, region) combination in parallel
  Int n_pairs    = length(vcf_pairs)
  Int n_snplists = length(MakeRegionSnplists.snplists)
  scatter (i in range(n_pairs * n_snplists)) {
    Int pair_idx    = i / n_snplists
    Int snplist_idx = i % n_snplists
    call FilterVCF {
      input:
        prefix  = vcf_pairs[pair_idx][0],
        vcf     = vcf_pairs[pair_idx][1],
        snplist = MakeRegionSnplists.snplists[snplist_idx]
    }
  }

  # Steps 3-6: per-VCF: concat -> subset query -> subset ref (parallel with FilterSNPs)
  scatter (i in range(n_pairs)) {
    Int first_chunk_idx = i * n_snplists
    call ConcatVCF {
      input:
        prefix   = vcf_pairs[i][0],
        all_vcfs = FilterVCF.filtered_vcf,
        disk_gb  = ceil(size(FilterVCF.filtered_vcf[first_chunk_idx], "GB") * n_snplists * 2) + 20
    }
    # Subset VCF to all shared HM3 variants
    call SubsetToPlink as SubsetQuery {
      input:
        input_files = [ConcatVCF.subset_vcf],
        snp_list    = plink_bim,
        out_prefix  = vcf_pairs[i][0] + "_query"
    }
    # Subset ref plink to same shared variants (uses query bim); runs parallel with FilterSNPs
    call SubsetToPlink as SubsetRef {
      input:
        input_files = [plink_bed, plink_bim, plink_fam],
        snp_list    = SubsetQuery.plink[1],
        out_prefix  = vcf_pairs[i][0] + "_ref"
    }
    # QC-filter query plink and build snplist; runs parallel with SubsetRef
    call FilterSNPs {
      input:
        plink_files = SubsetQuery.plink,
        prefix      = vcf_pairs[i][0],
        target_snps = target_snps
    }
    call PrepareDataset as PrepQuery {
      input:
        input_plink = SubsetQuery.plink,
        snp_list    = FilterSNPs.snplist,
        prefix      = vcf_pairs[i][0],
        max_het_F   = max_het_F,
        chunk_size  = chunk_size
    }
    call PrepareDataset as PrepRef {
      input:
        input_plink = SubsetRef.plink,
        snp_list    = FilterSNPs.snplist,
        prefix      = plink_prefix,
        max_het_F   = max_het_F,
        chunk_size  = chunk_size
    }
    call KingShards {
      input:
        query_chunks = PrepQuery.chunks,
        query_prefix = vcf_pairs[i][0],
        ref_chunks   = PrepRef.chunks,
        ref_prefix   = plink_prefix
    }
    call SummarizeKing {
      input:
        duplicates_con = KingShards.duplicates_con,
        query_fam      = PrepQuery.plink[2],
        query_prefix   = vcf_pairs[i][0],
        ref_prefix     = plink_prefix,
        het_excluded   = PrepQuery.excluded_samples
    }
  }

  call GatherResults {
    input:
      summaries  = SummarizeKing.summary,
      plots      = SummarizeKing.plot,
      out_prefix = out_prefix,
      aliases    = aliases
  }

  # =====================================================================
  # OPTIONAL: subset + rename to FinnGen IDs
  # =====================================================================
  if (run_rename) {
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
          resolved_mapping = GatherResults.id_mapping,
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
  }

  output {
    # ---- mandatory: kinship / ID resolution — always produced ----
    File combined_summary     = GatherResults.combined_summary
    File combined_plot        = GatherResults.combined_plot
    File id_mapping           = GatherResults.id_mapping
    File id_mapping_stats     = GatherResults.id_mapping_stats
    File id_mapping_md        = GatherResults.id_mapping_md
    File id_mapping_flowchart = GatherResults.id_mapping_flowchart

    # ---- optional: rename stage — only present when run_rename = true ----
    Array[File]? chrom_vcfs = ConcatChromVCF.out_vcf
    Array[File]? chrom_tbis = ConcatChromVCF.out_tbi
  }
}


# =========================================================================
# TASKS — kinship / duplicate detection (formerly exome_duplicates.wdl)
# =========================================================================

# -----------------------------------------------------------------------
# Step 1: Bin HM3 SNPs into Berisa LD blocks, greedily merge to N regions
# -----------------------------------------------------------------------
task MakeRegionSnplists {
  input {
    File    bim
    Int     n_regions     = 100
    File?   berisa_blocks
    String  berisa_url    = "https://bitbucket.org/nygcresearch/ldetect-data/raw/master/EUR/fourier_ls-all.bed"
  }

  command <<<
  set -euo pipefail

  BLOCKS="~{if defined(berisa_blocks) then select_first([berisa_blocks]) else ""}"
  BERISA_URL="~{berisa_url}"
  BIM="~{bim}"
  N_REGIONS=~{n_regions}

  if [[ -n "$BLOCKS" ]]; then
      cp "$BLOCKS" blocks.bed
  else
      curl -fsSL "$BERISA_URL" -o blocks.bed
  fi

  N_HM3=$(wc -l < "$BIM")
  N_BLOCKS=$(awk 'NR>1' blocks.bed | wc -l)
  echo "HM3 variants: $N_HM3  |  Berisa blocks: $N_BLOCKS  |  Target regions: ${N_REGIONS}"

  awk '
  NR==FNR {
      if (NR==1) next
      chr=$1; gsub(/^chr/,"",chr); start=$2+0; end=$3+0; idx=NR-2
      n = ++n_per_chr[chr]
      blk_start[chr,n]=start; blk_end[chr,n]=end; blk_idx[chr,n]=idx
      next
  }
  {
      chr=$1; gsub(/^chr/,"",chr); snp=$2; pos=$4+0
      n = n_per_chr[chr]
      if (!n) next
      lo=1; hi=n
      while (lo<=hi) {
          mid=int((lo+hi)/2)
          if      (blk_end[chr,mid]   < pos) lo=mid+1
          else if (blk_start[chr,mid] > pos) hi=mid-1
          else { print blk_idx[chr,mid], chr, pos, snp; break }
      }
  }
  ' blocks.bed "$BIM" | sort -k1,1n > assignments.txt

  N_ASSIGNED=$(wc -l < assignments.txt)
  echo "SNPs assigned: $N_ASSIGNED / $N_HM3"

  awk -v n_target="$N_REGIONS" -v total="$N_ASSIGNED" '
  BEGIN { regions_left=n_target; snps_left=total; target=int(total/n_target); region=0; cur=0; prev_block=-1 }
  {
      block=$1+0; chr=$2; pos=$3; snp=$4
      if (block != prev_block && cur >= target && region < n_target-1) {
          region++; snps_left -= cur; regions_left--
          target=int(snps_left/regions_left); cur=0
      }
      prev_block=block
      print chr "\t" pos "\t" snp >> (sprintf("snplist_%04d.txt", region))
      cur++
  }
  ' assignments.txt

  echo "Regions created: $(ls snplist_*.txt | wc -l)"
  >>>

  output {
    Array[File] snplists = glob("snplist_*.txt")
  }
}


# -----------------------------------------------------------------------
# Step 2: Extract one region from a VCF using bcftools -R / -T
# -----------------------------------------------------------------------
task FilterVCF {
  input {
    String vcf
    String prefix
    File   snplist
    Int    cpu       = 4
    Int    disk_gb   = 20
  }

  String snplist_id = basename(snplist, ".txt")
  String out_vcf    = prefix + "_" + snplist_id + ".vcf.gz"

  command <<<
  set -euo pipefail

  SNPLIST="~{snplist}"
  OUT_VCF="~{out_vcf}"
  VCF="~{vcf}"
  CPU=~{cpu}

  awk '{ print "chr" $1 "\t" $2 }' "$SNPLIST" | sort -k1,1V -k2,2n > positions.txt
  awk '$1 != prev { if (prev) print prev "\t" (lo-1) "\t" hi; prev=$1; lo=$2; hi=$2 }
       { if ($2<lo) lo=$2; if ($2>hi) hi=$2 }
       END { if (prev) print prev "\t" (lo-1) "\t" hi }' positions.txt > regions.bed
  export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
  bcftools view \
      --regions-file  regions.bed \
      --targets-file  positions.txt \
      --output-type z \
      --output        "$OUT_VCF" \
      --threads       "$CPU" \
      "$VCF"

  >>>

  output {
    File filtered_vcf = out_vcf
  }

  runtime {
    disks:  "local-disk ~{disk_gb} HDD"
    cpu:    cpu
  }
}


# -----------------------------------------------------------------------
# Step 3: Concatenate region chunks for one prefix, streaming from GCS
# -----------------------------------------------------------------------
task ConcatVCF {
  input {
    String        prefix
    Array[String] all_vcfs   # coerced from Array[File], not localized
    Int           disk_gb
    Int           cpu       = 4
    Int           memory_gb = 8
  }

  String out_vcf = prefix + "_subset.vcf.gz"

  command <<<
  set -euo pipefail

  PREFIX="~{prefix}"
  VCF_LIST_FILE="~{write_lines(all_vcfs)}"
  OUT_VCF="~{out_vcf}"
  CPU=~{cpu}

  export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

  grep -F "${PREFIX}_snplist_" "$VCF_LIST_FILE" \
      | awk -F'snplist_' '{split($2,a,"."); print a[1]+0, $0}' \
      | sort -k1,1n \
      | awk '{print $2}' > vcf_list.txt

  bcftools concat \
      --file-list vcf_list.txt \
      --output-type z \
      --output    "$OUT_VCF" \
      --threads   "$CPU"

  bcftools index -t --threads "$CPU" "$OUT_VCF"
  echo "Total variants: $(bcftools index -n "$OUT_VCF")"
  >>>

  output {
    File subset_vcf = out_vcf
    File subset_tbi = out_vcf + ".tbi"
  }

  runtime {
    memory: "~{memory_gb} GB"
    disks:  "local-disk ~{disk_gb} HDD"
    cpu:    cpu
  }
}


# -----------------------------------------------------------------------
# Step 4/6: Subset VCF or plink to a SNP list, return plink trio
# -----------------------------------------------------------------------
task SubsetToPlink {
  input {
    Array[String] input_files  # [vcf.gz] or [bed, bim, fam]
    File          snp_list
    String        out_prefix
    Int           cpu       = 16
    Int           disk_gb   = 50
  }

  Int disk_size = disk_gb
  Int memory_gb = 16

  command <<<
  set -euo pipefail

  SNP_LIST="~{snp_list}"
  OUT_PREFIX="~{out_prefix}"
  CPU=~{cpu}
  TOTAL_MEM_MB=12288

  to_fuse() { echo "$1" | sed 's|gs://[^/]*/|/mnt/disks/gcs/|'; }
  INPUT_FILE0=$(to_fuse "~{input_files[0]}")

  # normalise snp_list: bim files have 6 columns, extract col 2; plain ID files pass through
  awk 'NF > 1 {print $2} NF == 1 {print $1}' "$SNP_LIST" > _extract.txt

  if [[ "$INPUT_FILE0" == *.vcf.gz ]]; then
    INPUT_FLAG="--vcf $INPUT_FILE0 --double-id"
  else
    INPUT_FLAG="--bfile ${INPUT_FILE0%.bed}"
  fi

  plink2 \
    $INPUT_FLAG \
    --chr 1-22 \
    --extract _extract.txt \
    --rm-dup exclude-all \
    --make-bed \
    --out     "$OUT_PREFIX" \
    --threads "$CPU" \
    --memory  "$TOTAL_MEM_MB"

  echo "Variants: $(wc -l < "${OUT_PREFIX}.bim")"
  echo "Samples:  $(wc -l < "${OUT_PREFIX}.fam")"
  >>>

  output {
    Array[File] plink = [out_prefix + ".bed", out_prefix + ".bim", out_prefix + ".fam"]
  }

  runtime {
    memory: "~{memory_gb} GB"
    disks:  "local-disk ~{disk_size} HDD"
    cpu:    cpu
  }
}


# -----------------------------------------------------------------------
# Step 5: Select high-quality SNPs; pad with random SNPs to target count
# -----------------------------------------------------------------------
task FilterSNPs {
  input {
    Array[File] plink_files        # [bed, bim, fam] from SubsetToPlink
    String      prefix
    Int         target_snps
    Float       max_geno   = 0.05
    Float       hwe_p      = 0.0001
    Float       min_maf    = 0.01
    Int         cpu        = 4
    Int         memory_gb  = 8
    Int         disk_gb    = 20
  }

  command <<<
  set -euo pipefail

  BED="~{plink_files[0]}"
  OUT_PREFIX="~{prefix}"
  MAX_GENO=~{max_geno}
  HWE_P=~{hwe_p}
  MIN_MAF=~{min_maf}
  TARGET_SNPS=~{target_snps}

  PREFIX="${BED%.bed}"

  # Strict QC filter -> high-quality SNP list
  plink2 \
      --bfile      "$PREFIX" \
      --geno       "$MAX_GENO" \
      --hwe        "$HWE_P" 0 \
      --maf        "$MIN_MAF" \
      --write-snplist \
      --no-psam-pheno \
      --out        hq

  N_HQ=$(wc -l < hq.snplist)
  N_TOTAL=$(wc -l < "${PREFIX}.bim")
  echo "Total variants in plink:    $N_TOTAL"
  echo "High-quality SNPs after QC: $N_HQ  (geno<=${MAX_GENO}, HWE>${HWE_P}, MAF>=${MIN_MAF})"
  echo "Target SNP count:           ${TARGET_SNPS}"

  cp hq.snplist "${OUT_PREFIX}_hq_snplist.txt"

  # Append shuffled non-HQ variants below HQ, then take top target_snps
  awk '{print $2}' "${PREFIX}.bim" | sort > all_snps.txt
  sort hq.snplist | comm -23 all_snps.txt - | shuf >> hq.snplist
  head -n "$TARGET_SNPS" hq.snplist | sort -V > "${OUT_PREFIX}_snplist.txt"

  echo "Final SNP count: $(wc -l < "${OUT_PREFIX}_snplist.txt")"
  >>>

  output {
    File snplist    = "~{prefix}_snplist.txt"
    File hq_snplist = "~{prefix}_hq_snplist.txt"
  }

  runtime {
    memory: "~{memory_gb} GB"
    disks:  "local-disk ~{disk_gb} HDD"
    cpu:    cpu
  }
}


# -----------------------------------------------------------------------
# Steps 4+5: Subset to target SNPs, filter high-F samples, tag IDs, chunk
# -----------------------------------------------------------------------
task PrepareDataset {
  input {
    Array[File] input_plink
    File        snp_list
    String      prefix
    Float       max_het_F  = 0.3
    Int         chunk_size = 10000
    Int         cpu        = 16
    Int         memory_gb  = 8
  }

  String out_prefix = prefix + "_prepped"
  Int    disk_size  = ceil(size(input_plink[0], 'GB') * 4) + 20

  command <<<
  set -euo pipefail

  BED="~{input_plink[0]}"
  SNP_LIST="~{snp_list}"
  OUT_PREFIX="~{out_prefix}"
  PREFIX="~{prefix}"
  MAX_HET_F=~{max_het_F}
  CHUNK_SIZE=~{chunk_size}
  CPU=~{cpu}

  PLINK_PREFIX="${BED%.bed}"

  # 1. subset to target SNPs
  plink2 \
    --bfile   "$PLINK_PREFIX" \
    --extract "$SNP_LIST" \
    --make-bed \
    --out     subsetted \
    --threads "$CPU"

  # 2. het filter
  plink2 --bfile subsetted --freq --out subsetted --threads "$CPU"
  plink2 --bfile subsetted --read-freq subsetted.afreq --het --out subsetted --threads "$CPU"

  awk -v f="$MAX_HET_F" 'NR > 1 && $6 > f {print $1, $2}' subsetted.het > high_F_samples.txt
  echo "High-F samples excluded: $(wc -l < high_F_samples.txt)"

  awk -v f="$MAX_HET_F" 'NR == 1 || (NR > 1 && $6 > f)' subsetted.het > het_outliers.tsv

  plink2 \
    --bfile   subsetted \
    --remove  high_F_samples.txt \
    --make-bed \
    --out     "$OUT_PREFIX" \
    --threads "$CPU"

  # 3. annotate sample IDs
  awk -v p="${PREFIX}_" '{$2 = p $2; print}' "${OUT_PREFIX}.fam" > tmp.fam
  mv tmp.fam "${OUT_PREFIX}.fam"

  echo "Final samples:  $(wc -l < "${OUT_PREFIX}.fam")"
  echo "Final variants: $(wc -l < "${OUT_PREFIX}.bim")"

  # 4. split into sample chunks
  awk '{print $1, $2}' "${OUT_PREFIX}.fam" \
    | split -d -l "$CHUNK_SIZE" - "${OUT_PREFIX}_shard_"

  for shard in "${OUT_PREFIX}_shard_"*; do
    idx="${shard##*_shard_}"
    plink2 \
      --bfile   "$OUT_PREFIX" \
      --keep    "$shard" \
      --make-bed \
      --out     "${OUT_PREFIX}_chunk${idx}" \
      --threads "$CPU" \
      --silent
  done

  echo "Created $(ls "${OUT_PREFIX}_chunk"*.bed | wc -l) chunks of up to ${CHUNK_SIZE} samples"
  >>>

  output {
    Array[File] plink            = [out_prefix + ".bed", out_prefix + ".bim", out_prefix + ".fam"]
    File        excluded_samples = "het_outliers.tsv"
    Array[File] chunks           = glob(out_prefix + "_chunk*")
  }

  runtime {
    memory: "~{memory_gb} GB"
    disks:  "local-disk ~{disk_size} HDD"
    cpu:    cpu
  }
}


# -----------------------------------------------------------------------
# KING --duplicate across all N×M chunk pairs, sequential
# -----------------------------------------------------------------------
task KingShards {
  input {
    Array[File] query_chunks
    String      query_prefix
    Array[File] ref_chunks
    String      ref_prefix
    Int         cpu = 32
  }

  Int mem_raw   = ceil((size(query_chunks[0], 'GB') + size(ref_chunks[0], 'GB')) * 20)
  Int memory_gb = if mem_raw < 16 then 16 else mem_raw
  String out_prefix = query_prefix + "_vs_" + ref_prefix
  Int    disk_size  = ceil((size(query_chunks, 'GB') + size(ref_chunks, 'GB')) * 2) + 20

  command <<<
  set -euo pipefail

  QRY_PFX="~{query_prefix}_"
  QUERY_CHUNKS_FILE="~{write_lines(query_chunks)}"
  REF_CHUNKS_FILE="~{write_lines(ref_chunks)}"
  OUT_PREFIX="~{out_prefix}"

  grep '\.bed$' "$QUERY_CHUNKS_FILE" | sed 's/\.bed$//' > query_bases.txt
  grep '\.bed$' "$REF_CHUNKS_FILE"   | sed 's/\.bed$//' > ref_bases.txt

  echo "Query chunks: $(wc -l < query_bases.txt)"
  echo "Ref chunks:   $(wc -l < ref_bases.txt)"

  awk 'NR==FNR {a[$0]; next} {for (i in a) print i "," $0}' query_bases.txt ref_bases.txt > pairs.csv

  TOTAL_PAIRS=$(wc -l < pairs.csv)
  echo "Total chunk pairs: $TOTAL_PAIRS"

  JOB=0
  while IFS=',' read -r qbase rbase; do
    JOB=$((JOB + 1))
    echo "[$JOB/$TOTAL_PAIRS] $(basename "$qbase") vs $(basename "$rbase")"
    pfx="$(basename "$qbase")_vs_$(basename "$rbase")"
    king -b "${qbase}.bed,${rbase}.bed" --duplicate --cpu "$(nproc)" --prefix "$pfx"

    con="${pfx}.con"
    [[ ! -f "$con" ]] && continue

    awk -v qp="$QRY_PFX" -v hdr="$([ ! -f merged.con ] && echo 1 || echo 0)" '
      NR==1 { if (hdr) print; next }
      { if (($2 ~ "^"qp) != ($4 ~ "^"qp)) print }
    ' "$con" >> merged.con
  done < pairs.csv

  gzip -c merged.con > "${OUT_PREFIX}.con.gz"
  echo "Done. $(zcat "${OUT_PREFIX}.con.gz" | tail -n +2 | wc -l) duplicate pairs found."
  >>>

  output {
    File duplicates_con = out_prefix + ".con.gz"
  }

  runtime {
    memory: "~{memory_gb} GB"
    disks:  "local-disk ~{disk_size} HDD"
    cpu:    cpu
  }
}


# -----------------------------------------------------------------------
# Build per-sample summary from merged .con.gz
# -----------------------------------------------------------------------
task SummarizeKing {
  input {
    File   duplicates_con
    File   query_fam
    String query_prefix
    String ref_prefix
    File   het_excluded
  }

  Int    disk_size  = ceil(size(duplicates_con, 'GB') * 2) + 10
  String out_prefix = query_prefix + "_vs_" + ref_prefix

  command <<<
  set -euo pipefail

  QRY_PFX="~{query_prefix}_"
  OUT_PREFIX="~{out_prefix}"
  DUPLICATES_CON="~{duplicates_con}"
  QUERY_FAM="~{query_fam}"
  SUMMARY="${OUT_PREFIX}.summary.tsv"

  zcat "$DUPLICATES_CON" | awk -v qp="$QRY_PFX" '
    NR==1 { next }
    {
      if      ($2 ~ "^"qp) { qid=$1; rid=$3 }
      else if ($4 ~ "^"qp) { qid=$3; rid=$1 }
      else next
      m[qid] = (qid in m) ? m[qid] "," rid : rid
    }
    END { for (q in m) print q "\t" m[q] }
  ' > _matches.txt

  printf "QUERY\tDUPLICATES\n" > "$SUMMARY"
  awk '
    NR==FNR { m[$1]=$2; next }
    { print $1 "\t" ($1 in m ? m[$1] : "MISSING") }
  ' _matches.txt "$QUERY_FAM" >> "$SUMMARY"

  awk 'NR>1 {print $1 "\tHET_EXCLUDED"}' ~{het_excluded} >> "$SUMMARY"

  FOUND=$(awk 'NR>1 && $2!="MISSING" && $2!="HET_EXCLUDED"' "$SUMMARY" | wc -l)
  TOTAL=$(awk 'NR>1' "$SUMMARY" | wc -l)
  HET=$(awk 'NR>1 && $2=="HET_EXCLUDED"' "$SUMMARY" | wc -l)
  echo "$FOUND/$TOTAL query samples have duplicates in ref  ($HET excluded by het filter)"

  # Plot concordance diagnostics
  python3 << PYEOF
  import pandas as pd
  import matplotlib
  matplotlib.use("Agg")
  import matplotlib.pyplot as plt
  import matplotlib.gridspec as gridspec

  con_gz    = "$DUPLICATES_CON"
  out_png   = "${OUT_PREFIX}_concordance.png"
  threshold = 0.9

  df      = pd.read_csv(con_gz, sep="\t", compression="gzip")
  n_pairs = len(df)
  n_dup   = (df["Concord"] >= threshold).sum()
  is_dup  = df["Concord"] >= threshold
  colors  = ["#d62728" if d else "#aec7e8" for d in is_dup]

  fig = plt.figure(figsize=(14, 5))
  fig.suptitle(f"{con_gz}  |  {n_pairs} pairs  |  {n_dup} duplicates (Concord >= {threshold})",fontsize=10, y=1.01)
  gs = gridspec.GridSpec(1, 3, wspace=0.35)

  ax1 = fig.add_subplot(gs[0])
  ax1.hist(df["Concord"], bins=100, color="#4878d0", edgecolor="none")
  ax1.axvline(threshold, color="#d62728", linestyle="--", linewidth=1)
  ax1.set_xlabel("Concordance"); ax1.set_ylabel("Pairs"); ax1.set_title("Concordance distribution")

  ax2 = fig.add_subplot(gs[1])
  ax2.scatter(df["N_IBS0"], df["Concord"], c=colors, s=6, alpha=0.6, linewidths=0)
  ax2.set_xlabel("N_IBS0"); ax2.set_ylabel("Concordance"); ax2.set_title("IBS0 vs Concordance")

  ax3 = fig.add_subplot(gs[2])
  ax3.hist(df["N"], bins=50, color="#4878d0", edgecolor="none")
  ax3.set_xlabel("N SNPs"); ax3.set_ylabel("Pairs"); ax3.set_title("SNP count per pair")

  fig.tight_layout()
  fig.savefig(out_png, dpi=150, bbox_inches="tight")
  print(f"Saved: {out_png}  ({n_dup}/{n_pairs} pairs above threshold)")
  PYEOF
  >>>

  output {
    File summary = out_prefix + ".summary.tsv"
    File plot    = out_prefix + "_concordance.png"
  }

  runtime {
    disks: "local-disk ~{disk_size} HDD"
  }
}


# -----------------------------------------------------------------------
# Gather: combine all summaries and plots into one (always runs)
# -----------------------------------------------------------------------
task GatherResults {
  input {
    Array[File] summaries
    Array[File] plots
    String      out_prefix
    File?       aliases
    String      docker = "eu.gcr.io/finngen-refinery-dev/exome_bioinf:dup.2"
  }

  command <<<
  set -euo pipefail
  PREFIX="~{out_prefix}"

  # ── 1. Combine per-dataset summaries into combined_summary.tsv ──────────
  mapfile -t summary_files < "~{write_lines(summaries)}"
  # Write header once with DATASET column prepended
  awk 'NR==1{print "DATASET\t" $0; exit}' "${summary_files[0]}" > "${PREFIX}_combined_summary.tsv"
  # Append data rows from each dataset
  for sf in "${summary_files[@]}"; do
      dataset=$(basename "$sf" | sed 's/\.summary\.tsv$//' | sed 's/_vs_.*//')
      awk -v ds="$dataset" 'NR>1{print ds "\t" $0}' "$sf"
  done >> "${PREFIX}_combined_summary.tsv"

  # ── 2. Stack per-dataset concordance plots ──────────────────────────────
  python3 -c "import matplotlib,os; matplotlib.use('Agg'); import matplotlib.pyplot as P,matplotlib.image as I; pf=[l.strip() for l in open('~{write_lines(plots)}') if l.strip()]; fig,ax=P.subplots(len(pf),1,figsize=(14,5*len(pf))); ax=[ax] if len(pf)==1 else list(ax); [a.imshow(I.imread(f)) or a.axis('off') or a.set_title(os.path.basename(f).replace('_concordance.png',''),fontsize=12) for a,f in zip(ax,pf)]; P.tight_layout(); P.savefig('${PREFIX}_concordance.png',dpi=150,bbox_inches='tight')"

  # ── 3. Resolve mapping, stats, and flowchart ────────────────────────────
  python3 /scripts/resolve_mapping.py \
    "${PREFIX}_combined_summary.tsv" \
    ~{if defined(aliases) then "--aliases " + select_first([aliases]) else ""} \
    --out "${PREFIX}_id_mapping.tsv"

  >>>

  output {
    File combined_summary = out_prefix + "_combined_summary.tsv"
    File combined_plot    = out_prefix + "_concordance.png"
    File id_mapping       = out_prefix + "_id_mapping.tsv"
    File id_mapping_stats = out_prefix + "_id_mapping_stats.tsv"
    File id_mapping_md    = out_prefix + "_id_mapping_stats.md"
    File id_mapping_flowchart = out_prefix + "_id_mapping_flowchart.png"
  }

  meta {
    volatile: true
  }

  runtime {
    docker: docker
    disks:  "local-disk 20 HDD"
  }
}


# =========================================================================
# TASKS — subset + rename to FinnGen IDs (formerly exome_rename.wdl)
# =========================================================================

# -----------------------------------------------------------------------
# Query CHROM+POS for all datasets for one chrom in parallel (cached).
# Output: one {dataset}_{chrom}_pos.txt per dataset.
# -----------------------------------------------------------------------
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


# -----------------------------------------------------------------------
# Build size-based chunk regions for all datasets × chroms.
# Queries file stats (fast — index only), then loops over position files.
# Output: Array[Array[String]] tasks with columns vcf | prefix | region
# prefix = DATASET_CHUNKIDX e.g. botnia_001
# -----------------------------------------------------------------------
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


# -----------------------------------------------------------------------
# Stream one region, subset + rename samples.
# prefix = DATASET_CHUNKIDX → dataset derived by cutting at '_'
# -----------------------------------------------------------------------
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

  export HTS_HTTP_VERSION=1.1
  cut -f1 mapping.tsv > sample_list.txt

  success=0
  for attempt in 1 2 3; do
    export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
    rm -f subset.vcf.gz
    if bcftools view -r "~{region}" -t "~{region}" --samples-file sample_list.txt --force-samples --threads ~{cpu} -Oz -o subset.vcf.gz "~{vcf}"; then
      success=1
      break
    fi
    echo "Attempt $attempt failed for ~{region}" >&2
    [[ $attempt -lt 3 ]] && sleep 15
  done
  [[ $success -eq 1 ]] || { echo "All 3 attempts failed for ~{region}" >&2; exit 1; }

  # Reheader locally — no GCS, no token risk
  bcftools reheader --samples mapping.tsv --output "~{out}" subset.vcf.gz
  rm subset.vcf.gz

  echo "DONE ~{out} samples=$(bcftools query -l ~{out} | wc -l)"
  >>>

  output {
    File out_vcf = out
  }

  runtime {
    cpu:   cpu
    disks: "local-disk ~{disk_gb} HDD"
  }
}


# -----------------------------------------------------------------------
# Concat all chunks for one dataset+chrom into a single VCF.
# Receives all chunk paths as strings (no localization), greps to select
# matching chunks, sorts by chunk index, concatenates.
# -----------------------------------------------------------------------
task ConcatChromVCF {
  input {
    Array[String] all_vcf_paths
    String        dataset
    String        vcf
    String        chrom
    String        suffix
    Int           chrom_size_mb
    Int           cpu     = 16
    Int           disk_gb = chrom_size_mb * 2 / 1024 + 10
  }

  String base = sub(basename(vcf, ".vcf.gz"), "\\.QC_ANNOTATED$", "")
  String out  = base + ".QC_ANNOTATED_" + suffix + "_" + chrom + ".vcf.gz"

  command <<<
  set -euo pipefail

  grep "~{dataset}_~{chrom}_" ~{write_lines(all_vcf_paths)} | \
    awk -F/ '{print $NF "\t" $0}' | sort -k1,1 | cut -f2- > chunks.txt
  echo "~{dataset} ~{chrom}: $(wc -l < chunks.txt) chunks" >&2
  export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
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
