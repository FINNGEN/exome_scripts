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

    # Step 4+5: subset to target SNPs, het-filter, annotate IDs
    call PrepareDataset as PrepQuery {
      input:
        input_plink = SubsetQuery.plink,
        snp_list    = SharedSNPs.snp_list,
        prefix      = pair[0],
        max_het_F   = max_het_F
    }

    call PrepareDataset as PrepRef {
      input:
        input_plink = SubsetRef.plink,
        snp_list    = SharedSNPs.snp_list,
        prefix      = plink_prefix,
        max_het_F   = max_het_F
    }

    call KingShards {
      input:
        query_plink  = PrepQuery.plink,
        query_prefix = pair[0],
        ref_plink    = PrepRef.plink,
        ref_prefix   = plink_prefix,
        chunk_size   = chunk_size
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
# Step 3: Downsample shared SNP list to target (0 = use all)
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
# Steps 1+2+4+5: Subset VCF or plink to SNP list, return plink (no annotation)
# -----------------------------------------------------------------------
task SubsetToPlink {
  input {
    String      input_type   # "vcf" or "plink"
    Array[File] input_files  # [vcf] or [bed, bim, fam]
    File        snp_list
    String      out_prefix
    Int         cpu       = 8
    Int         memory_gb = 16
  }

  Int disk_size = ceil(size(input_files[0], 'GB') * 3) + 20

  command <<<
  set -euo pipefail
  TOTAL_MEM_MB=$(awk '/MemAvailable/ {print int($2/1024)}' /proc/meminfo)

  # normalise snp_list: bim files have 6 columns, extract col 2; plain ID files pass through
  awk 'NF > 1 {print $2} NF == 1 {print $1}' "~{snp_list}" > _extract.txt

  if [[ "~{input_type}" == "vcf" ]]; then
    INPUT_FLAG="--vcf ~{input_files[0]}"
  else
    INPUT_FLAG="--bfile ~{sub(input_files[0], '\\.bed$', '')}"
  fi

  plink2 \
    $INPUT_FLAG \
    --chr 1-22 \
    --extract _extract.txt \
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
# Steps 4+5: Subset to target SNPs, filter high-F samples, tag sample IDs
# -----------------------------------------------------------------------
task PrepareDataset {
  input {
    Array[File] input_plink
    File        snp_list
    String      prefix
    Float       max_het_F = 0.3
    Int         cpu       = 4
    Int         memory_gb = 16
  }

  String out_prefix = prefix + "_prepped"
  Int    disk_size  = ceil(size(input_plink[0], 'GB') * 2) + 20

  command <<<
  set -euo pipefail

  # 1. subset to target SNPs
  plink2 \
    --bed "~{input_plink[0]}" --bim "~{input_plink[1]}" --fam "~{input_plink[2]}" \
    --extract "~{snp_list}" \
    --make-bed \
    --out subsetted \
    --threads ~{cpu}

  # 2. het filter on subsetted data
  plink2 --bfile subsetted --freq --out subsetted --threads ~{cpu}
  plink2 --bfile subsetted --read-freq subsetted.afreq --het --out subsetted --threads ~{cpu}

  awk 'NR > 1 && $6 > ~{max_het_F} {print $1, $2}' subsetted.het > high_F_samples.txt
  echo "High-F samples excluded: $(wc -l < high_F_samples.txt)"

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
  >>>

  output {
    Array[File] plink            = [out_prefix + ".bed", out_prefix + ".bim", out_prefix + ".fam"]
    File        excluded_samples = "high_F_samples.txt"
  }

  runtime {
    memory: "~{memory_gb} GB"
    disks:  "local-disk ~{disk_size} HDD"
    cpu:    cpu
  }
}


# -----------------------------------------------------------------------
# Step 8: Sharded KING --duplicate, output merged gzipped .con
# -----------------------------------------------------------------------
task KingShards {
  input {
    Array[File] query_plink
    String      query_prefix
    Array[File] ref_plink
    String      ref_prefix
    Int         chunk_size
    Int         cpu       = 32
    Int         memory_gb = 64
  }

  String out_prefix = query_prefix + "_vs_" + ref_prefix
  Int    disk_size  = ceil((size(query_plink[0], 'GB') + size(ref_plink[0], 'GB')) * 3) + 20

  command <<<
  set -euo pipefail

  QUERY_PLINK="~{sub(query_plink[0], '\\.bed$', '')}"
  REF_PLINK="~{sub(ref_plink[0], '\\.bed$', '')}"
  QRY_PFX="~{query_prefix}_"
  CHUNK_SIZE=~{chunk_size}

  echo "Query: $QUERY_PLINK  ($(wc -l < "~{query_plink[2]}") samples)"
  echo "Ref:   $REF_PLINK  ($(wc -l < "~{ref_plink[2]}") samples)"
  echo "Chunk size: $CHUNK_SIZE"
  echo ""

  TMP_DIR="/dev/shm/king_~{out_prefix}"
  mkdir -p "$TMP_DIR"
  trap "rm -rf '$TMP_DIR'" EXIT

  awk '{print $1, $2}' "${REF_PLINK}.fam" \
    | split -d -l "$CHUNK_SIZE" - "${TMP_DIR}/chunk_"
  mapfile -t CHUNKS < <(ls "${TMP_DIR}/chunk_"*)
  echo "Running KING across ${#CHUNKS[@]} shards..."

  PIPELINE_SH="${TMP_DIR}/pipeline.sh"
  WORKDIR=$(pwd)
  for chunk in "${CHUNKS[@]}"; do
    echo "plink2 --bfile ${WORKDIR}/${REF_PLINK} --keep ${chunk} --make-bed --out ${chunk}_ref --threads 1 --memory 4000 --silent \
      && king -b ${WORKDIR}/${QUERY_PLINK}.bed,${chunk}_ref.bed \
              --bim ${WORKDIR}/${QUERY_PLINK}.bim,${chunk}_ref.bim \
              --fam ${WORKDIR}/${QUERY_PLINK}.fam,${chunk}_ref.fam \
              --duplicate --cpu 1 --prefix ${chunk}_king > ${chunk}.kinglog 2>&1 \
      && rm -f ${chunk}_ref.bed ${chunk}_ref.bim ${chunk}_ref.fam ${chunk}_ref.log \
      || echo FAILED > ${chunk}.failed"
  done > "$PIPELINE_SH"
  parallel --bar -j "$(nproc)" < "$PIPELINE_SH"

  mapfile -t FAILED < <(ls "${TMP_DIR}/chunk_"*.failed 2>/dev/null || true)
  if [[ ${#FAILED[@]} -gt 0 ]]; then
    echo "WARNING: ${#FAILED[@]} shard(s) failed:"
    for f in "${FAILED[@]}"; do echo "  $(basename "$f" .failed)"; done
  fi

  HEADER_WRITTEN=0
  for chunk in "${CHUNKS[@]}"; do
    con="${chunk}_king.con"
    [[ ! -f "$con" ]] && continue
    if [[ $HEADER_WRITTEN -eq 0 ]]; then
      head -1 "$con" > merged.con
      HEADER_WRITTEN=1
    fi
    awk -v qp="$QRY_PFX" 'NR > 1 {
      qry1 = ($2 ~ "^"qp); qry2 = ($4 ~ "^"qp)
      if (qry1 != qry2) print
    }' "$con" >> merged.con
  done
  [[ $HEADER_WRITTEN -eq 0 ]] && printf "FID1\tID1\tFID2\tID2\tN_SNP\tKinship\tIBS0\n" > merged.con

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
    Int    memory_gb = 4
  }

  Int    disk_size  = ceil(size(duplicates_con, 'GB') * 2) + 10
  String out_prefix = query_prefix + "_vs_" + ref_prefix

  command <<<
  set -euo pipefail

  QRY_PFX="~{query_prefix}_"
  REF_PFX="~{ref_prefix}_"
  SUMMARY="~{out_prefix}.summary.tsv"

  # pass 1: build qid -> comma-separated rids from con.gz
  zcat "~{duplicates_con}" | awk -v qp="$QRY_PFX" -v rp="$REF_PFX" '
    NR==1 { next }
    {
      id1=$2; id2=$4
      if (substr(id1,1,length(qp))==qp) {
        qid=substr(id1,length(qp)+1)
        rid=(substr(id2,1,length(rp))==rp) ? substr(id2,length(rp)+1) : id2
      } else if (substr(id2,1,length(qp))==qp) {
        qid=substr(id2,length(qp)+1)
        rid=(substr(id1,1,length(rp))==rp) ? substr(id1,length(rp)+1) : id1
      } else next
      m[qid] = (qid in m) ? m[qid] "," rid : rid
    }
    END { for (q in m) print q "\t" m[q] }
  ' > _matches.txt

  # pass 2: join against fam order, write summary
  printf "QUERY\tDUPLICATES\n" > "$SUMMARY"
  awk -v qp="$QRY_PFX" '
    NR==FNR { m[$1]=$2; next }
    { raw=$2; qid=(substr(raw,1,length(qp))==qp) ? substr(raw,length(qp)+1) : raw
      print qid "\t" (qid in m ? m[qid] : "MISSING") }
  ' _matches.txt "~{query_fam}" >> "$SUMMARY"

  FOUND=$(awk 'NR>1 && $2!="MISSING"' "$SUMMARY" | wc -l)
  TOTAL=$(awk 'NR>1' "$SUMMARY" | wc -l)
  echo "$FOUND/$TOTAL query samples have duplicates in ref"
  >>>

  output {
    File summary = out_prefix + ".summary.tsv"
  }

  runtime {
    memory: "~{memory_gb} GB"
    disks:  "local-disk ~{disk_size} HDD"
  }
}
