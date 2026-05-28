version 1.0

workflow exome_king {
  input {
    Array[Array[String]] vcf_pairs
    File   plink_bed
    String plink_prefix
    Int    chunk_size  = 10000
    Int    target_snps
    Float  max_het_F   = 0.3
  }

  File plink_bim = sub(plink_bed, "\\.bed$", ".bim")
  File plink_fam = sub(plink_bed, "\\.bed$", ".fam")

  scatter (pair in vcf_pairs) {

    # Step 1+2: build cached full-shared-variant datasets
    call SubsetToPlink as SubsetQuery {
      input:
        input_type  = "vcf",
        input_files = [pair[1]],
        snp_list    = plink_bim,
        out_prefix  = pair[0] + "_shared"
    }

    call SubsetToPlink as SubsetRef {
      input:
        input_type  = "plink",
        input_files = [plink_bed, plink_bim, plink_fam],
        snp_list    = SubsetQuery.plink[1],
        out_prefix  = plink_prefix + "_shared"
    }

    # Step 3: downsample to target SNPs
    call SharedSNPs {
      input:
        bim         = SubsetQuery.plink[1],
        target_snps = target_snps
    }

    # Step 4+5: subset to target SNPs, het-filter, annotate IDs, chunk
    call PrepareDataset as PrepQuery {
      input:
        input_plink = SubsetQuery.plink,
        snp_list    = SharedSNPs.snp_list,
        prefix      = pair[0],
        max_het_F   = max_het_F,
        chunk_size  = chunk_size
    }

    call PrepareDataset as PrepRef {
      input:
        input_plink = SubsetRef.plink,
        snp_list    = SharedSNPs.snp_list,
        prefix      = plink_prefix,
        max_het_F   = max_het_F,
        chunk_size  = chunk_size
    }

    call KingShards {
      input:
        query_chunks = PrepQuery.chunks,
        query_prefix = pair[0],
        ref_chunks   = PrepRef.chunks,
        ref_prefix   = plink_prefix
    }

    call SummarizeKing {
      input:
        duplicates_con = KingShards.duplicates_con,
        query_fam      = PrepQuery.plink[2],
        query_prefix   = pair[0],
        ref_prefix     = plink_prefix
    }
  }

  output {
    Array[File] summary                = SummarizeKing.summary
    Array[File] duplicates_con         = KingShards.duplicates_con
    Array[File] excluded_samples_query = PrepQuery.excluded_samples
    Array[File] excluded_samples_ref   = PrepRef.excluded_samples
  }
}


# -----------------------------------------------------------------------
# Step 3: Downsample shared SNP list to target
# -----------------------------------------------------------------------
task SharedSNPs {
  input {
    File bim
    Int  target_snps
    Int  memory_gb = 2
  }

  command <<<
  set -euo pipefail

  awk '{print $2}' "~{bim}" > all_snps.txt
  echo "Shared variants: $(wc -l < all_snps.txt)"
  shuf -n ~{target_snps} all_snps.txt > snp_list.txt
  echo "SNP list size: $(wc -l < snp_list.txt)"
  >>>

  output {
    File snp_list = "snp_list.txt"
  }

  runtime {
    memory: "~{memory_gb} GB"
    disks:  "local-disk 10 HDD"
  }
}


# -----------------------------------------------------------------------
# Steps 1+2: Subset VCF or plink to SNP list, return plink (no annotation)
# -----------------------------------------------------------------------
task SubsetToPlink {
  input {
    String      input_type   # "vcf" or "plink"
    Array[File] input_files  # [vcf] or [bed, bim, fam]
    File        snp_list
    String      out_prefix
    Int         cpu       = 16
  }

  Int disk_size = ceil(size(input_files[0], 'GB') * 3) + 20
  Int memory_gb = 16

  command <<<
  set -euo pipefail
  # Cap plink2 at 12 GB — it will spill to disk rather than OOM the VM
  TOTAL_MEM_MB=12288

  # normalise snp_list: bim files have 6 columns, extract col 2; plain ID files pass through
  awk 'NF > 1 {print $2} NF == 1 {print $1}' "~{snp_list}" > _extract.txt

  if [[ "~{input_type}" == "vcf" ]]; then
    INPUT_FLAG="--vcf ~{input_files[0]} --double-id"
  else
    INPUT_FLAG="--bfile ~{sub(input_files[0], '\\.bed$', '')}"
  fi

  plink2 \
    $INPUT_FLAG \
    --chr 1-22 \
    --extract _extract.txt \
    --rm-dup exclude-all \
    --make-bed \
    --out "~{out_prefix}" \
    --threads ~{cpu} \
    --memory "$TOTAL_MEM_MB"

  echo "Variants: $(wc -l < "~{out_prefix}.bim")"
  echo "Samples:  $(wc -l < "~{out_prefix}.fam")"
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

  # 1. subset to target SNPs
  plink2 \
    --bfile "~{sub(input_plink[0], '\\.bed$', '')}" \
    --extract "~{snp_list}" \
    --make-bed \
    --out subsetted \
    --threads ~{cpu}

  # 2. het filter
  plink2 --bfile subsetted --freq --out subsetted --threads ~{cpu}
  plink2 --bfile subsetted --read-freq subsetted.afreq --het --out subsetted --threads ~{cpu}

  awk 'NR > 1 && $6 > ~{max_het_F} {print $1, $2}' subsetted.het > high_F_samples.txt
  echo "High-F samples excluded: $(wc -l < high_F_samples.txt)"

  # Write het stats (header + outlier rows) for downstream inspection
  awk 'NR == 1 || (NR > 1 && $6 > ~{max_het_F})' subsetted.het > het_outliers.tsv

  plink2 \
    --bfile subsetted \
    --remove high_F_samples.txt \
    --make-bed \
    --out "~{out_prefix}" \
    --threads ~{cpu}

  # 3. annotate sample IDs last
  awk -v p="~{prefix}_" '{$2 = p $2; print}' "~{out_prefix}.fam" > tmp.fam
  mv tmp.fam "~{out_prefix}.fam"

  echo "Final samples:  $(wc -l < "~{out_prefix}.fam")"
  echo "Final variants: $(wc -l < "~{out_prefix}.bim")"

  # 4. split into sample chunks
  awk '{print $1, $2}' "~{out_prefix}.fam" \
    | split -d -l ~{chunk_size} - "~{out_prefix}_shard_"

  for shard in "~{out_prefix}_shard_"*; do
    idx="${shard##*_shard_}"
    plink2 \
      --bfile "~{out_prefix}" \
      --keep "$shard" \
      --make-bed \
      --out "~{out_prefix}_chunk${idx}" \
      --threads ~{cpu} \
      --silent
  done

  echo "Created $(ls "~{out_prefix}_chunk"*.bed | wc -l) chunks of up to ~{chunk_size} samples"
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
# Step 8: KING --duplicate across all N×M chunk pairs, sequential
# -----------------------------------------------------------------------
task KingShards {
  input {
    Array[File] query_chunks   # all chunk bed+bim+fam files for query
    String      query_prefix
    Array[File] ref_chunks     # all chunk bed+bim+fam files for ref
    String      ref_prefix
    Int         cpu = 32
  }

  Int mem_raw  = ceil((size(query_chunks[0], 'GB') + size(ref_chunks[0], 'GB')) * 20)
  Int memory_gb = if mem_raw < 16 then 16 else mem_raw
  String out_prefix = query_prefix + "_vs_" + ref_prefix
  Int    disk_size  = ceil((size(query_chunks, 'GB') + size(ref_chunks, 'GB')) * 2) + 20

  command <<<
  set -euo pipefail

  QRY_PFX="~{query_prefix}_"

  # Write chunk file lists; filter to .bed and strip extension to get basepaths
  grep '\.bed$' "~{write_lines(query_chunks)}" | sed 's/\.bed$//' > query_bases.txt
  grep '\.bed$' "~{write_lines(ref_chunks)}"   | sed 's/\.bed$//' > ref_bases.txt

  echo "Query chunks: $(wc -l < query_bases.txt)"
  echo "Ref chunks:   $(wc -l < ref_bases.txt)"

  # Generate N×M pairs CSV (query_base,ref_base)
  awk 'NR==FNR {a[$0]; next} {for (i in a) print i "," $0}' query_bases.txt ref_bases.txt > pairs.csv

  TOTAL_PAIRS=$(wc -l < pairs.csv)
  echo "Total chunk pairs: $TOTAL_PAIRS"

  # Run KING on each pair sequentially with all CPUs
  HEADER_WRITTEN=0
  JOB=0
  while IFS=',' read -r qbase rbase; do
    JOB=$((JOB + 1))
    echo "[$JOB/$TOTAL_PAIRS] $(basename "$qbase") vs $(basename "$rbase")"
    pfx="$(basename "$qbase")_vs_$(basename "$rbase")"
    king -b "${qbase}.bed,${rbase}.bed" --duplicate --cpu "$(nproc)" --prefix "$pfx"

    con="${pfx}.con"
    [[ ! -f "$con" ]] && continue

    cat "$con" | awk -v qp="$QRY_PFX" -v hdr="$([ ! -f merged.con ] && echo 1 || echo 0)" '
      NR==1 { if (hdr) print; next }
      { if (($2 ~ "^"qp) != ($4 ~ "^"qp)) print }
    ' >> merged.con
  done < pairs.csv

  gzip -c merged.con > "~{out_prefix}.con.gz"
  echo "Done. $(zcat "~{out_prefix}.con.gz" | tail -n +2 | wc -l) duplicate pairs found."
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
# Step 9: Build per-sample summary from merged .con.gz
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
  SUMMARY="~{out_prefix}.summary.tsv"

  # pass 1: FID ($1/$3) is the original sample name; use IID ($2/$4) only to identify query vs ref
  zcat "~{duplicates_con}" | awk -v qp="$QRY_PFX" '
    NR==1 { next }
    {
      if      ($2 ~ "^"qp) { qid=$1; rid=$3 }
      else if ($4 ~ "^"qp) { qid=$3; rid=$1 }
      else next
      m[qid] = (qid in m) ? m[qid] "," rid : rid
    }
    END { for (q in m) print q "\t" m[q] }
  ' > _matches.txt

  # pass 2: join against fam using FID ($1) directly — no prefix stripping needed
  printf "QUERY\tDUPLICATES\n" > "$SUMMARY"
  awk '
    NR==FNR { m[$1]=$2; next }
    { print $1 "\t" ($1 in m ? m[$1] : "MISSING") }
  ' _matches.txt "~{query_fam}" >> "$SUMMARY"

  FOUND=$(awk 'NR>1 && $2!="MISSING"' "$SUMMARY" | wc -l)
  TOTAL=$(awk 'NR>1' "$SUMMARY" | wc -l)
  echo "$FOUND/$TOTAL query samples have duplicates in ref"
  >>>

  output {
    File summary = out_prefix + ".summary.tsv"
  }

  runtime {
    disks:  "local-disk ~{disk_size} HDD"
  }
}
