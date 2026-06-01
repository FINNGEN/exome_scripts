version 1.0

workflow exome_duplicates {
  input {
    Array[Array[String]] vcf_pairs
    File   plink_bed
    String plink_prefix
    Int    n_regions     = 100
    File?  berisa_blocks
    File?  aliases
    Int    target_snps   = 10000
    Float  max_het_F     = 0.3
    Int    chunk_size    = 10000
  }

  File plink_bim = sub(plink_bed, "\\.bed$", ".bim")
  File plink_fam = sub(plink_bed, "\\.bed$", ".fam")

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
        ref_prefix     = plink_prefix
    }
  }

  call GatherResults {
    input:
      summaries    = SummarizeKing.summary,
      plots        = SummarizeKing.plot,
      snplists     = FilterSNPs.snplist,
      plink_prefix = plink_prefix,
      aliases      = aliases
  }

  output {
    Array[File]        subset_vcfs            = ConcatVCF.subset_vcf
    Array[File]        subset_tbis            = ConcatVCF.subset_tbi
    Array[File]        snplists               = FilterSNPs.snplist
    Array[File]        hq_snplists            = FilterSNPs.hq_snplist
    Array[File]        summary                = SummarizeKing.summary
    Array[File]        concordance_plots      = SummarizeKing.plot
    Array[File]        duplicates_con         = KingShards.duplicates_con
    Array[File]        excluded_samples_query = PrepQuery.excluded_samples
    Array[File]        excluded_samples_ref   = PrepRef.excluded_samples
    File               combined_plot          = GatherResults.combined_plot
    File               resolved_mapping       = GatherResults.resolved_mapping
    File               resolved_stats_tsv     = GatherResults.resolved_stats_tsv
    File               resolved_stats_md      = GatherResults.resolved_stats_md

  }
}


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

  export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

  awk '{ print "chr" $1 "\t" $2 }' "$SNPLIST" | sort -k1,1V -k2,2n > positions.txt

  awk '$1 != prev { if (prev) print prev "\t" (lo-1) "\t" hi; prev=$1; lo=$2; hi=$2 }
       { if ($2<lo) lo=$2; if ($2>hi) hi=$2 }
       END { if (prev) print prev "\t" (lo-1) "\t" hi }' positions.txt > regions.bed

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
    Array[File] input_files  # [vcf.gz] or [bed, bim, fam]
    File        snp_list
    String      out_prefix
    Int         cpu       = 16
  }

  Int disk_size = ceil(size(input_files[0], 'GB') * 3) + 20
  Int memory_gb = 16

  command <<<
  set -euo pipefail

  SNP_LIST="~{snp_list}"
  INPUT_FILE0="~{input_files[0]}"
  OUT_PREFIX="~{out_prefix}"
  CPU=~{cpu}
  TOTAL_MEM_MB=12288

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

  FOUND=$(awk 'NR>1 && $2!="MISSING"' "$SUMMARY" | wc -l)
  TOTAL=$(awk 'NR>1' "$SUMMARY" | wc -l)
  echo "$FOUND/$TOTAL query samples have duplicates in ref"

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
    Array[File] snplists
    String      plink_prefix
    File?       aliases
  }

  command <<<
  set -euo pipefail

  SUMMARIES_FILE="~{write_lines(summaries)}"
  PLOTS_FILE="~{write_lines(plots)}"
  SNPLISTS_FILE="~{write_lines(snplists)}"

  python3 << PYEOF
import pandas as pd
import os

with open("$SUMMARIES_FILE") as f:
    summary_files = [l.strip() for l in f if l.strip()]
with open("$SNPLISTS_FILE") as f:
    snplist_files = [l.strip() for l in f if l.strip()]

dfs = []
global_rows = []

for sf, sl in zip(summary_files, snplist_files):
    dataset = os.path.basename(sf).replace(".summary.tsv", "").split("_vs_")[0]
    df = pd.read_csv(sf, sep="\t")
    df.insert(0, "DATASET", dataset)
    dfs.append(df)

    n_total     = len(df)
    n_matched   = (df["DUPLICATES"] != "MISSING").sum()
    n_ambiguous = df["DUPLICATES"].apply(lambda x: isinstance(x, str) and "," in x).sum()
    n_snps      = sum(1 for _ in open(sl))

    global_rows.append({
        "DATASET":     dataset,
        "TOTAL_QUERY": n_total,
        "N_MATCHED":   n_matched,
        "PCT_MATCHED": f"{n_matched / n_total * 100:.1f}%",
        "N_AMBIGUOUS": n_ambiguous,
        "N_SNPS":      n_snps,
    })

combined = pd.concat(dfs, ignore_index=True)
combined.to_csv("~{plink_prefix}_EXOME_summary.tsv", sep="\t", index=False)

global_df = pd.DataFrame(global_rows)
global_df.to_csv("~{plink_prefix}_EXOME_global_summary.tsv", sep="\t", index=False)

print(global_df.to_string(index=False))
PYEOF

  python3 << PYEOF
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os

with open("$PLOTS_FILE") as f:
    plot_files = [l.strip() for l in f if l.strip()]

fig, axes = plt.subplots(len(plot_files), 1, figsize=(14, 5 * len(plot_files)))
if len(plot_files) == 1:
    axes = [axes]
for ax, pf in zip(axes, plot_files):
    img = mpimg.imread(pf)
    ax.imshow(img)
    ax.axis("off")
    ax.set_title(os.path.basename(pf).replace("_concordance.png", ""), fontsize=12, pad=8)
plt.tight_layout()
plt.savefig("~{plink_prefix}_EXOME_concordance.png", dpi=150, bbox_inches="tight")
print(f"Combined {len(plot_files)} plots")
PYEOF

  python3 << PYEOF
import os, random
import pandas as pd

PREFIX     = "~{plink_prefix}_EXOME"
ALIAS_PATH = "~{if defined(aliases) then select_first([aliases]) else ""}"

RESOLVED           = {"ID_CONFIRMED","RESOLVED_BY_ID","RESOLVED_BY_ALIAS","INFERRED_BY_ELIMINATION","UNIQUE"}
CONFLICT_PRIORITY  = {"ID_CONFIRMED":0,"RESOLVED_BY_ALIAS":1,"RESOLVED_BY_ID":2,"INFERRED_BY_ELIMINATION":3,"UNIQUE":4}
SUMMARY_GROUPS = [
    ("MATCHED", "samples with a final QRY→REF mapping in the output", [
        ("ID_CONFIRMED",            "single candidate; KING match confirmed by matching IDs"),
        ("RESOLVED_BY_ID",          "twins in ref; query ID matched one candidate"),
        ("RESOLVED_BY_ALIAS",       "twins in ref; candidates are known aliases of each other"),
        ("INFERRED_BY_ELIMINATION", "twins in ref; all other candidates already claimed"),
        ("UNIQUE",                  "single candidate; matched by genetics only"),
        ("CONFLICT_KEPT",           "contested ref ID; kept after priority tiebreak"),
    ]),
    ("DROPPED", "found by KING but excluded from final mapping", [
        ("CONFLICT_DROPPED",     "contested ref ID; lost tiebreak; REF_MAPPED = NA"),
        ("AMBIGUOUS_UNRESOLVED", "multiple ref candidates; no resolution possible"),
        ("AMBIGUOUS_ALL_TAKEN",  "multiple ref candidates; all already claimed"),
    ]),
    ("NO MATCH", "absent from ref or below KING concordance threshold", [
        ("MISSING", "no KING match found"),
    ]),
]

def load_aliases(path):
    if not path or not os.path.exists(path): return {}
    groups, n = {}, 0
    with open(path) as f:
        for line in f:
            ids = [x.strip() for x in line.strip().split("\t") if x.strip()]
            if len(ids) < 2: continue
            g = frozenset(ids)
            for id_ in ids: groups[id_] = (g, ids[0])
            n += 1
    print(f"Loaded {n} alias groups ({len(groups)} IDs)")
    return groups

def _alias_resolve(query, cands, ag):
    qi = ag.get(query); ai = [ag.get(c) for c in cands]; nk = sum(1 for i in ai if i is not None)
    if qi is not None:
        m = [c for c,x in zip(cands,ai) if x is not None and x[0]==qi[0]]
        if len(m)==1: return m[0], ""
    if nk==len(cands)>0 and len({x[0] for x in ai})==1:
        canon = ai[0][1]; return (canon if canon in cands else cands[0]), ""
    note = "no_candidates_in_alias_file" if nk==0 else \
           f"{nk}/{len(cands)}_candidates_in_alias_file" if nk<len(cands) else \
           "candidates_in_different_alias_groups"
    return None, note

def _row(ds, q, ref, cands, status, note=""):
    return dict(DATASET=ds, QUERY=q, REF_MAPPED=ref, CANDIDATES=cands, STATUS=status, ALIAS_NOTE=note)

def initial_categorise(df, ag):
    records = []
    for _, row in df.iterrows():
        ds, query, raw = row["DATASET"], str(row["QUERY"]).strip(), str(row["DUPLICATES"]).strip()
        if raw in ("MISSING","nan",""): records.append(_row(ds,query,"NA","","MISSING")); continue
        seen, cands = set(), []
        for c in raw.split(","):
            c=c.strip()
            if c and c not in seen: cands.append(c); seen.add(c)
        if len(cands)==1:
            ref=cands[0]
            if query==ref: records.append(_row(ds,query,ref,raw,"ID_CONFIRMED"))
            elif ag:
                resolved,_ = _alias_resolve(query,cands,ag)
                records.append(_row(ds,query,ref,raw,"RESOLVED_BY_ALIAS" if resolved else "UNIQUE"))
            else: records.append(_row(ds,query,ref,raw,"UNIQUE"))
        else:
            if query in cands: records.append(_row(ds,query,query,raw,"RESOLVED_BY_ID"))
            elif ag:
                resolved,note = _alias_resolve(query,cands,ag)
                records.append(_row(ds,query,resolved or "AMBIGUOUS",raw,"RESOLVED_BY_ALIAS" if resolved else "AMBIGUOUS_UNRESOLVED",note))
            else: records.append(_row(ds,query,"AMBIGUOUS",raw,"AMBIGUOUS_UNRESOLVED"))
    return pd.DataFrame(records)

def disambiguate_by_elimination(result):
    claimed = set(result.loc[result["STATUS"].isin(RESOLVED),"REF_MAPPED"])
    changed = True
    while changed:
        changed = False
        for idx, row in result[result["STATUS"]=="AMBIGUOUS_UNRESOLVED"].iterrows():
            free = [c.strip() for c in row["CANDIDATES"].split(",") if c.strip() and c.strip() not in claimed]
            if len(free)==1: result.at[idx,"REF_MAPPED"]=free[0]; result.at[idx,"STATUS"]="INFERRED_BY_ELIMINATION"; claimed.add(free[0]); changed=True
            elif len(free)==0: result.at[idx,"STATUS"]="AMBIGUOUS_ALL_TAKEN"; changed=True
    return result

def check_surjectivity(result, rng):
    ref_to_idx = {}
    for idx, row in result[result["STATUS"].isin(RESOLVED)].iterrows():
        ref_to_idx.setdefault(row["REF_MAPPED"],[]).append(idx)
    for _, indices in ref_to_idx.items():
        if len(indices)==1: continue
        def pri(i): return CONFLICT_PRIORITY.get(result.at[i,"STATUS"].split("[")[0],99)
        best = min(pri(i) for i in indices)
        winner = rng.choice([i for i in indices if pri(i)==best])
        for idx in indices:
            orig = result.at[idx,"STATUS"]
            if idx==winner: result.at[idx,"STATUS"]=f"CONFLICT_KEPT[{orig}]"
            else: result.at[idx,"STATUS"]=f"CONFLICT_DROPPED[{orig}]"; result.at[idx,"REF_MAPPED"]="NA"
    return result

def _pct(n,t): return f"{n/t*100:.1f}%" if t else "n/a"

def build_stats_table(result):
    n, datasets = len(result), sorted(result["DATASET"].unique())
    gc = result["STATUS"].value_counts()
    dsc = {ds: result[result["DATASET"]==ds]["STATUS"].value_counts() for ds in datasets}
    def pfx(p,vc): return sum(v for k,v in vc.items() if k.startswith(p))
    cols = ["SECTION","GROUP","STATUS","TOTAL"]+datasets+["PCT","NOTES"]
    def blank(): return {c:"" for c in cols}
    tot, brk = [], []
    for gh, gn, statuses in SUMMARY_GROUPS:
        prefixes = [s for s,_ in statuses]; gcount = sum(pfx(p,gc) for p in prefixes)
        row = dict(SECTION="TOTALS",GROUP=gh,STATUS="",TOTAL=gcount,PCT=_pct(gcount,n),NOTES=gn)
        for ds in datasets: row[ds]=sum(pfx(p,dsc[ds]) for p in prefixes)
        tot.append(row)
        for prefix, desc in statuses:
            cnt = pfx(prefix,gc)
            if cnt==0: continue
            if prefix=="CONFLICT_KEPT":
                nc = result[result["STATUS"].str.startswith("CONFLICT_KEPT")]["REF_MAPPED"].nunique()
                desc += f"; {nc} ref IDs contested, avg {gcount/nc:.1f} queries/ref" if nc else ""
            row = dict(SECTION="BREAKDOWN",GROUP=gh,STATUS=prefix,TOTAL=cnt,PCT=_pct(cnt,n),NOTES=desc)
            for ds in datasets: row[ds]=pfx(prefix,dsc[ds])
            brk.append(row)
    grand = dict(SECTION="TOTALS",GROUP="TOTAL",STATUS="",TOTAL=n,PCT="100.0%",NOTES="")
    for ds in datasets: grand[ds]=(result["DATASET"]==ds).sum()
    tot.append(grand)
    sep=blank(); sep["SECTION"]="---"
    return pd.DataFrame(tot+[sep]+brk,columns=cols)

def df_to_md(df):
    cols=list(df.columns)
    rows=["| "+" | ".join(str(c) for c in cols)+" |","| "+" | ".join("---" for _ in cols)+" |"]
    for _,row in df.iterrows(): rows.append("| "+" | ".join("" if str(v)=="nan" else str(v) for v in row)+" |")
    return "\n".join(rows)

ag     = load_aliases(ALIAS_PATH)
df     = pd.read_csv(f"{PREFIX}_summary.tsv", sep="\t")
result = initial_categorise(df, ag)
result = disambiguate_by_elimination(result)
result = check_surjectivity(result, random.Random(42))

result[["QUERY","REF_MAPPED","DATASET","STATUS","CANDIDATES","ALIAS_NOTE"]].to_csv(f"{PREFIX}_resolved.tsv", sep="\t", index=False)

stats     = build_stats_table(result)
stats.to_csv(f"{PREFIX}_resolved_stats.tsv", sep="\t", index=False)

totals    = stats[stats["SECTION"]=="TOTALS"].drop(columns=["SECTION","STATUS"]).reset_index(drop=True)
breakdown = stats[stats["SECTION"]=="BREAKDOWN"].drop(columns="SECTION").reset_index(drop=True)
with open(f"{PREFIX}_resolved_stats.md","w") as fh:
    fh.write(f"## Mapping Totals\n\n{df_to_md(totals)}\n\n## Mapping Breakdown\n\n{df_to_md(breakdown)}\n")

n=len(result); gc=result["STATUS"].value_counts()
def pfx(p,vc): return sum(v for k,v in vc.items() if k.startswith(p))
print(f"\nResolved mapping ({PREFIX}): {n:,} samples")
for gh,_,statuses in SUMMARY_GROUPS:
    gn=sum(pfx(s,gc) for s,_ in statuses)
    print(f"  {gh:<10} {gn:>7,}  ({_pct(gn,n)})")
PYEOF
  >>>

  output {
    File combined_plot       = plink_prefix + "_EXOME_concordance.png"
    File resolved_mapping    = plink_prefix + "_EXOME_resolved.tsv"
    File resolved_stats_tsv  = plink_prefix + "_EXOME_resolved_stats.tsv"
    File resolved_stats_md   = plink_prefix + "_EXOME_resolved_stats.md"
  }

  meta {
    volatile: true
  }

  runtime {
    disks: "local-disk 20 HDD"
  }
}
