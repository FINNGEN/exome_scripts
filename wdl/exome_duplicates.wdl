version 1.0

workflow exome_duplicates {
  input {
    Array[Array[String]] vcf_pairs
    File   plink_bed
    String plink_prefix
    File?  bim
    Int    chunk_size = 10000
  }

  File plink_bim = sub(plink_bed, "\\.bed$", ".bim")
  Array[File] plink_input_files = [plink_bed, plink_bim, sub(plink_bed, "\\.bed$", ".fam")]

  scatter (pair in vcf_pairs) {
    # pre-filter VCF to bim SNPs (parallel bcftools per chrom, streaming)
    call SubsetVCF {
      input:
        prefix    = pair[0],
        input_vcf = pair[1],
        bim       = select_first([bim, plink_bim])
    }

    # subset plink to snplist and rename sample IDs with plink_prefix
    call PlinkSubset {
      input:
        dataset_prefix = pair[0],
        plink_prefix   = plink_prefix,
        plink_files    = plink_input_files,
        snplist        = SubsetVCF.snplist
    }

    # find cross-dataset duplicates by genotype concordance
    call RunGtcheck {
      input:
        query_vcf   = SubsetVCF.filtered_vcf,
        query_tbi   = SubsetVCF.filtered_vcf_tbi,
        plink_files = PlinkSubset.plink_out,
        prefix      = pair[0] + "_gtcheck",
        chunk_size  = chunk_size
    }
  }

  output {
    Array[File] gtcheck_raw     = RunGtcheck.gtcheck_raw
    Array[File] gtcheck_summary = RunGtcheck.gtcheck_summary
  }
}

task SubsetVCF {
  input {
    String prefix
    File   input_vcf
    File   bim
    Int    cpu       = 24
    Int    memory_gb = cpu
  }

  File   input_vcf_tbi = input_vcf + ".tbi"
  String output_vcf    = prefix + ".subset.vcf.gz"
  Int    disk_size     = ceil(size(input_vcf, 'GB') * 3)

  command <<<
  set -euo pipefail
  VCF="~{input_vcf}"
  touch "~{input_vcf_tbi}"
  BIM="~{bim}"
  OUTPUT_VCF="~{output_vcf}"
  CHUNKS=$(( $(nproc) - 1 ))
  if [[ $CHUNKS -lt 1 ]]; then CHUNKS=1; fi

  # chromosomes present in the bim, normalised to chr-prefix
  mkdir -p ./tmp
  mapfile -t chromosomes < <(awk '{print $1}' "$BIM" | sort -u | awk '{print (/^chr/ ? $0 : "chr"$0)}')

  N_BEFORE=$(bcftools index -s "$VCF" | awk '{sum+=$3} END {print sum}')
  echo "Variants before filter: $N_BEFORE"
  echo "Processing ${#chromosomes[@]} chromosomes in parallel ($CHUNKS jobs)..."
  echo ""

  for chrom in "${chromosomes[@]}"; do
    safe_chrom=$(echo "$chrom" | sed 's/[*:\/]/_/g')
    awk -v c="$chrom" '$1 == c || "chr"$1 == c {print c"\t"$4}' "$BIM" > "./tmp/pos_${safe_chrom}.txt"
  done

  # build sample rename map: OLD_NAME -> PREFIX_OLD_NAME
  bcftools query -l "$VCF" | awk -v p="~{prefix}" '{print $0"\t"p"_"$0}' > sample_rename.txt

  SCRIPT_DIR=$(mktemp -d)
  for chrom in "${chromosomes[@]}"; do
    safe_chrom=$(echo "$chrom" | sed 's/[*:\/]/_/g')
    cat > "${SCRIPT_DIR}/run_${safe_chrom}.sh" << SCRIPT
#!/bin/bash
bcftools view -r "${chrom}" -T "./tmp/pos_${safe_chrom}.txt" "${VCF}" | bcftools reheader -s ./sample_rename.txt | bgzip -c > "chunk_${safe_chrom}.vcf.gz"
echo "Done: ${chrom}"
SCRIPT
  done

  ls "${SCRIPT_DIR}"/run_*.sh | parallel -j "$CHUNKS" 'bash {}'

  for chrom in "${chromosomes[@]}"; do
    safe_chrom=$(echo "$chrom" | sed 's/[*:\/]/_/g')
    echo "chunk_${safe_chrom}.vcf.gz"
  done > chunk_list.txt

  bcftools concat -n -Oz -o "$OUTPUT_VCF" -f chunk_list.txt
  rm -f chunk_*.vcf.gz chunk_list.txt

  bcftools index -t "$OUTPUT_VCF"
  N_AFTER=$(bcftools index -s "$OUTPUT_VCF" | awk '{sum+=$3} END {print sum}')
  echo "Variants after filter: $N_AFTER"

  bcftools query -f '%ID\n' "$OUTPUT_VCF" > "~{prefix}.snplist.txt"
  >>>

  output {
    File filtered_vcf     = output_vcf
    File filtered_vcf_tbi = output_vcf + ".tbi"
    File snplist          = prefix + ".snplist.txt"
  }

  runtime {
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_size} HDD"
    cpu: cpu
  }
}

task PlinkSubset {
  input {
    String      dataset_prefix
    String      plink_prefix
    Array[File] plink_files
    File        snplist
    Int         cpu       = 4
    Int         memory_gb = 16
  }

  File   plink_bed  = plink_files[0]
  Int    disk_size  = ceil(size(plink_bed, 'GB') * 2) + 10
  String out_prefix = dataset_prefix + "_" + plink_prefix

  command <<<
  set -euo pipefail
  PLINK_PREFIX="~{sub(plink_bed, '\\.bed$', '')}"

  echo "Ref:      $PLINK_PREFIX"
  echo "Variants: $(wc -l < '~{snplist}') SNPs"
  echo ""

  # build renamed fam before extract: old_FID old_IID new_FID new_IID
  awk -v p="~{plink_prefix}" 'BEGIN{OFS="\t"} {print $1,$2,$1,p"_"$2}' "$PLINK_PREFIX.fam" > update_ids.txt
  plink2 \
    --bfile "$PLINK_PREFIX" \
    --update-ids update_ids.txt \
    --make-just-fam \
    --out renamed_tmp \
    --threads ~{cpu}

  plink2 \
    --bfile "$PLINK_PREFIX" \
    --fam renamed_tmp.fam \
    --extract "~{snplist}" \
    --make-bed \
    --out "~{out_prefix}" \
    --threads ~{cpu}

  echo "Done."
  echo "  Samples:  $(wc -l < '~{out_prefix}.fam')"
  echo "  Variants: $(wc -l < '~{out_prefix}.bim')"
  >>>

  output {
    Array[File] plink_out = ["~{out_prefix}.bed", "~{out_prefix}.bim", "~{out_prefix}.fam"]
  }

  runtime {
    memory: "~{memory_gb} GB"
    disks:  "local-disk ~{disk_size} HDD"
    cpu:    cpu
  }
}

task RunGtcheck {
  input {
    File        query_vcf
    File        query_tbi
    Array[File] plink_files
    String      prefix
    Int         chunk_size = 100
    Int         cpu        = 32
    Int         memory_gb  = 32
  }

  File plink_bed = plink_files[0]
  Int  disk_size = ceil(size(query_vcf, 'GB') + size(plink_bed, 'GB')) * 3 + 20

  command <<<
  set -euo pipefail
  touch "~{query_tbi}"

  # symlink plink files to a consistent local prefix
  ln -s "~{plink_files[0]}" plink_in.bed
  ln -s "~{plink_files[1]}" plink_in.bim
  ln -s "~{plink_files[2]}" plink_in.fam

  QUERY_VCF="~{query_vcf}"
  PREFIX="~{prefix}"

  echo "Query: $QUERY_VCF ($(bcftools query -l "$QUERY_VCF" | wc -l) samples)"
  echo "Ref:   $(wc -l < plink_in.fam) samples, $(wc -l < plink_in.bim) variants"
  echo "Chunks: ~{chunk_size} samples each"
  echo ""

  # fam already has renamed IIDs from PlinkSubset — split directly
  awk '{print $1, $2}' plink_in.fam | split -d -l ~{chunk_size} - "${PREFIX}_chunk_"
  mapfile -t CHUNKS < <(ls "${PREFIX}_chunk_"*)

  NJOBS=$(( ${#CHUNKS[@]} < $(nproc) ? ${#CHUNKS[@]} : $(nproc) ))
  MEM_PER_JOB=$(( ~{memory_gb} * 1024 / NJOBS ))

  # split.sh: plink2 conversion per chunk
  SPLIT_SH="${PREFIX}_split.sh"
  for chunk in "${CHUNKS[@]}"; do
    echo "plink2 --bfile plink_in --keep ${chunk} --export vcf id-paste=iid bgz --output-chr chrM --out ${chunk}_ref --threads 1 --memory ${MEM_PER_JOB}"
  done > "$SPLIT_SH"
  parallel -j "$(nproc)" < "$SPLIT_SH"

  # index chunk_00 and copy tbi to all other chunks
  bcftools index -t "${PREFIX}_chunk_00_ref.vcf.gz"
  for chunk in "${CHUNKS[@]}"; do
    [[ "$chunk" == "${PREFIX}_chunk_00" ]] && continue
    cp "${PREFIX}_chunk_00_ref.vcf.gz.tbi" "${chunk}_ref.vcf.gz.tbi"
  done

  # gtcheck.sh: one gtcheck per chunk
  GTCHECK_SH="${PREFIX}_gtcheck.sh"
  for chunk in "${CHUNKS[@]}"; do
    echo "bcftools gtcheck --no-HWE-prob -g ${chunk}_ref.vcf.gz ${QUERY_VCF} > ${chunk}.gtcheck && rm ${chunk}_ref.vcf.gz ${chunk}_ref.vcf.gz.tbi"
  done > "$GTCHECK_SH"
  parallel -j "$(nproc)" < "$GTCHECK_SH"

  cat "${PREFIX}_chunk_"*.gtcheck > "${PREFIX}.gtcheck"
  rm -f "${PREFIX}_chunk_"* "$SPLIT_SH" "$GTCHECK_SH"

  python3 - "${PREFIX}.gtcheck" << 'PYEOF' > "${PREFIX}.summary.tsv"
import sys

stats = {}
with open(sys.argv[1]) as f:
    for line in f:
        if not line.startswith('DCv2'):
            continue
        parts = line.rstrip('\n').split('\t')
        query, ref = parts[1], parts[2]
        nsites = int(parts[5])
        rate = float(parts[3]) / nsites if nsites > 0 else 1.0
        if query not in stats:
            stats[query] = {'best': (rate, ref), 'second': None, 'sum': rate, 'count': 1}
        else:
            s = stats[query]
            s['sum'] += rate
            s['count'] += 1
            if rate < s['best'][0]:
                s['second'] = s['best']
                s['best'] = (rate, ref)
            elif s['second'] is None or rate < s['second'][0]:
                s['second'] = (rate, ref)

print('\t'.join(['QUERY','BEST_MATCH','BEST_RATE','2ND_MATCH','2ND_RATE','AVG_OTHERS','RATIO']))
rows = []
for query, s in stats.items():
    best_rate, best_ref = s['best']
    second_rate, second_ref = s['second'] if s['second'] else (float('nan'), 'N/A')
    avg_others = (s['sum'] - best_rate) / (s['count'] - 1) if s['count'] > 1 else float('nan')
    ratio = best_rate / avg_others if avg_others > 0 else float('nan')
    rows.append((best_rate, query, best_ref, second_ref, second_rate, avg_others, ratio))
rows.sort()
for best_rate, query, best_ref, second_ref, second_rate, avg_others, ratio in rows:
    print(f'{query}\t{best_ref}\t{best_rate:.6f}\t{second_ref}\t{second_rate:.6f}\t{avg_others:.6f}\t{ratio:.6f}')
PYEOF

  echo ""
  echo "Done."
  echo "  Sites compared: $(awk '/^DCv2/ {print $6; exit}' "${PREFIX}.gtcheck")"
  echo "  Likely duplicates (ratio < 0.1): $(awk -F'\t' 'NR>1 && $7 < 0.1' "${PREFIX}.summary.tsv" | wc -l)"
  >>>

  output {
    File gtcheck_raw     = prefix + ".gtcheck"
    File gtcheck_summary = prefix + ".summary.tsv"
  }

  runtime {
    memory: "~{memory_gb} GB"
    disks:  "local-disk ~{disk_size} HDD"
    cpu:    cpu
  }
}
