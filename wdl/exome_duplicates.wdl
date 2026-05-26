version 1.0

workflow exome_duplicates {
  input {
    Array[Array[String]] vcf_pairs
    File   plink_bed
    String plink_prefix
    File?  bim
    Int    chunk_size   = 10000
    Int    target_snps
  }

  File plink_bim = sub(plink_bed, "\\.bed$", ".bim")
  Array[File] plink_input_files = [plink_bed, plink_bim, sub(plink_bed, "\\.bed$", ".fam")]

  scatter (pair in vcf_pairs) {
    # pre-filter VCF to bim SNPs (parallel bcftools per chrom, streaming)
    call SubsetVCF {
      input:
        prefix      = pair[0],
        input_vcf   = pair[1],
        bim         = select_first([bim, plink_bim]),
        target_snps = target_snps
    }

    # subset plink to snplist and rename sample IDs with plink_prefix
    call PlinkSubset {
      input:
        dataset_prefix = pair[0],
        plink_prefix   = plink_prefix,
        plink_files    = plink_input_files,
        snplist        = SubsetVCF.snplist
    }

    # split ref VCF into indexed chunks (cached independently)
    call SplitRefVCF {
      input:
        ref_vcf    = PlinkSubset.ref_vcf,
        ref_tbi    = PlinkSubset.ref_tbi,
        chunk_size = chunk_size
    }

    # find cross-dataset duplicates by genotype concordance
    call RunGtcheck {
      input:
        query_vcf = SubsetVCF.filtered_vcf,
        query_tbi = SubsetVCF.filtered_vcf_tbi,
        ref_vcfs  = SplitRefVCF.chunk_vcfs,
        ref_tbis  = SplitRefVCF.chunk_tbis,
        prefix    = pair[0] + "_gtcheck"
    }

    call SummarizeGtcheck {
      input:
        gtcheck_raw = RunGtcheck.gtcheck_raw,
        prefix      = pair[0] + "_gtcheck"
    }
  }

  output {
    Array[File] gtcheck_raw     = RunGtcheck.gtcheck_raw
    Array[File] gtcheck_summary = SummarizeGtcheck.gtcheck_summary
    Array[File] gtcheck_plot    = SummarizeGtcheck.gtcheck_plot
  }
}

task SubsetVCF {
  input {
    String prefix
    File   input_vcf
    File   bim
    Int    target_snps
    Int    cpu       = 24
    Int    memory_gb = cpu
  }

  File   input_vcf_tbi = input_vcf + ".tbi"
  String output_vcf    = prefix + ".subset.vcf.gz"
  Int    disk_size     = ceil(size(input_vcf, 'GB') * 3)

  command <<<
  set -euo pipefail
  VCF="~{input_vcf}"
  INPUT_VCF_TBI="~{input_vcf_tbi}"
  BIM="~{bim}"
  OUTPUT_VCF="~{output_vcf}"
  PREFIX="~{prefix}"
  touch "$INPUT_VCF_TBI"
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
  bcftools query -l "$VCF" | awk -v p="$PREFIX" '{print $0"\t"p"_"$0}' > sample_rename.txt

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

  bcftools query -f '%ID\n' "$OUTPUT_VCF" > "${PREFIX}.snplist.txt"
  TARGET_SNPS=~{target_snps}

  # intersect VCF IDs with BIM IDs to ensure ID consistency after position-based filtering
  awk '{print $2}' "$BIM" > bim_ids.txt
  comm -12 <(sort "${PREFIX}.snplist.txt") <(sort bim_ids.txt) > "${PREFIX}.snplist.intersect.txt"
  N_INTERSECT=$(wc -l < "${PREFIX}.snplist.intersect.txt")
  echo "Variants in BIM ID intersection: $N_INTERSECT"

  if [[ $N_INTERSECT -gt $TARGET_SNPS ]]; then
    shuf -n "$TARGET_SNPS" "${PREFIX}.snplist.intersect.txt" | sort -V > "${PREFIX}.snplist.txt"
  else
    sort -V "${PREFIX}.snplist.intersect.txt" > "${PREFIX}.snplist.txt"
  fi

  bcftools view -i "ID=@${PREFIX}.snplist.txt" -Oz -o "${PREFIX}.subset2.vcf.gz" "$OUTPUT_VCF"
  mv "${PREFIX}.subset2.vcf.gz" "$OUTPUT_VCF"
  bcftools index -t "$OUTPUT_VCF"
  echo "SNPs after intersection+subsample: $(wc -l < "${PREFIX}.snplist.txt")"
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
  Int    disk_size  = ceil(size(plink_bed, 'GB') * 3) + 10
  String out_prefix = dataset_prefix + "_" + plink_prefix

  command <<<
  set -euo pipefail
  PLINK_PREFIX="~{sub(plink_bed, '\\.bed$', '')}"
  SNPLIST="~{snplist}"
  RENAME_PREFIX="~{plink_prefix}"
  OUT_PREFIX="~{out_prefix}"
  CPU=~{cpu}

  echo "Ref:      $PLINK_PREFIX"
  echo "Variants: $(wc -l < "$SNPLIST") SNPs"
  echo ""

  # build renamed fam before extract: old_FID old_IID new_FID new_IID
  awk -v p="$RENAME_PREFIX" 'BEGIN{OFS="\t"} {print $1,$2,$1,p"_"$2}' "$PLINK_PREFIX.fam" > update_ids.txt
  plink2 \
    --bfile "$PLINK_PREFIX" \
    --update-ids update_ids.txt \
    --make-just-fam \
    --out renamed_tmp \
    --threads $CPU

  plink2 \
    --bfile "$PLINK_PREFIX" \
    --fam renamed_tmp.fam \
    --extract "$SNPLIST" \
    --make-bed \
    --out "$OUT_PREFIX" \
    --threads $CPU

  plink2 \
    --bfile "$OUT_PREFIX" \
    --export vcf id-paste=iid bgz \
    --output-chr chrM \
    --out "$OUT_PREFIX" \
    --threads $CPU
  bcftools index -t "${OUT_PREFIX}.vcf.gz"

  echo "Done."
  echo "  Samples:  $(wc -l < "${OUT_PREFIX}.fam")"
  echo "  Variants: $(wc -l < "${OUT_PREFIX}.bim")"
  >>>

  output {
    Array[File] plink_out = ["~{out_prefix}.bed", "~{out_prefix}.bim", "~{out_prefix}.fam"]
    File        ref_vcf   = out_prefix + ".vcf.gz"
    File        ref_tbi   = out_prefix + ".vcf.gz.tbi"
  }

  runtime {
    memory: "~{memory_gb} GB"
    disks:  "local-disk ~{disk_size} HDD"
    cpu:    cpu
  }
}

task SplitRefVCF {
  input {
    File   ref_vcf
    File   ref_tbi
    Int    chunk_size
    Int    cpu        = 8
    Int    memory_gb  = cpu*2
  }

  Int disk_size = ceil(size(ref_vcf, 'GB') * 2) + 10

  command <<<
  set -euo pipefail
  REF_VCF="~{ref_vcf}"
  REF_TBI="~{ref_tbi}"
  CHUNK_SIZE=~{chunk_size}
  touch "$REF_TBI"

  N_SAMPLES=$(bcftools query -l "$REF_VCF" | wc -l)
  N_CHUNKS=$(( (N_SAMPLES + CHUNK_SIZE - 1) / CHUNK_SIZE ))
  echo "Splitting $N_SAMPLES samples into $N_CHUNKS chunks of up to $CHUNK_SIZE..."

  bcftools query -l "$REF_VCF" | split -d -l $CHUNK_SIZE - "samples_chunk_"
  mapfile -t slist_files < <(ls samples_chunk_*)

  INPUT_SIZE=$(stat -c%s "$REF_VCF")
  (
    while true; do
      sleep 10
      DONE=$(stat -c%s samples_chunk_*.bcf 2>/dev/null | awk '{s+=$1} END {print s+0}')
      PCT=$(( INPUT_SIZE > 0 ? DONE * 100 / INPUT_SIZE : 0 ))
      echo "  written: $(( DONE / 1024 / 1024 )) / $(( INPUT_SIZE / 1024 / 1024 )) MB  (${PCT}%)"
    done
  ) &
  PROGRESS_PID=$!

  for chunk in samples_chunk_*; do
    echo "bcftools view -S ${chunk} -Ob --write-index -o ${chunk}.bcf ${REF_VCF}"
  done | parallel -j "$(nproc)"

  kill $PROGRESS_PID 2>/dev/null
  wait $PROGRESS_PID 2>/dev/null || true
  rm -f "${slist_files[@]}"
  echo "Done: $N_CHUNKS BCF chunks created."
  >>>

  output {
    Array[File] chunk_vcfs = glob("samples_chunk_*.bcf")
    Array[File] chunk_tbis = glob("samples_chunk_*.bcf.csi")
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
    Array[File] ref_vcfs
    Array[File] ref_tbis
    String      prefix
    Int         cpu = 32
  }

  Int effective_cpu = if length(ref_vcfs) < cpu then length(ref_vcfs) else cpu
  Int memory_gb     =  ceil(effective_cpu/2.0)
  Int disk_size     = ceil(size(query_vcf, 'GB') + size(ref_vcfs[0], 'GB') * length(ref_vcfs)) + 20

  command <<<
  set -euo pipefail
  QUERY_VCF="~{query_vcf}"
  QUERY_TBI="~{query_tbi}"
  PREFIX="~{prefix}"
  REF_VCFS_LIST="~{write_lines(ref_vcfs)}"
  REF_TBIS_LIST="~{write_lines(ref_tbis)}"
  touch "$QUERY_TBI"

  echo "Query: $QUERY_VCF ($(bcftools query -l "$QUERY_VCF" | wc -l) samples)"
  echo "Chunks: $(wc -l < "$REF_VCFS_LIST") VCF chunks"
  echo ""

  PIPELINE_SH="${PREFIX}_pipeline.sh"
  paste "$REF_VCFS_LIST" "$REF_TBIS_LIST" | while IFS=$'\t' read -r vcf tbi; do
    mv "$tbi" "${vcf}.csi"
    name=$(basename "$vcf" .bcf)
    echo "bcftools gtcheck --no-HWE-prob -g ${vcf} ${QUERY_VCF} > chunk_${name}.gtcheck"
  done > "$PIPELINE_SH"
  parallel -j "$(nproc)" < "$PIPELINE_SH"

  cat chunk_*.gtcheck > "${PREFIX}.gtcheck"
  rm -f chunk_*.gtcheck "$PIPELINE_SH"

  echo "Done. $(awk '/^DCv2/' "${PREFIX}.gtcheck" | wc -l) pairwise comparisons written."
  >>>

  output {
    File gtcheck_raw = prefix + ".gtcheck"
  }

  runtime {
    memory: "~{memory_gb} GB"
    disks:  "local-disk ~{disk_size} HDD"
    cpu:    effective_cpu
  }
}

task SummarizeGtcheck {
  input {
    File   gtcheck_raw
    String prefix
  }


  command <<<
  set -euo pipefail
  GTCHECK_RAW="~{gtcheck_raw}"
  PREFIX="~{prefix}"

  python3 - "$GTCHECK_RAW" "${PREFIX}.summary.tsv" "${PREFIX}.distribution.png" "$PREFIX" << 'PYEOF'
import sys, math, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

gtcheck_file, summary_file, out_png, prefix = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

stats = {}
with open(gtcheck_file) as f:
    for line in f:
        if not line.startswith('DCv2'):
            continue
        parts = line.rstrip('\n').split('\t')
        query, ref = parts[1], parts[2]
        nsites = int(parts[5])
        rate = float(parts[3]) / nsites if nsites > 0 else 1.0
        if query not in stats:
            near = [(rate, ref)] if rate < 0.02 else []
            stats[query] = {'best': (rate, ref), 'second': None, 'sum': rate, 'count': 1, 'near': near}
        else:
            s = stats[query]
            s['sum'] += rate
            s['count'] += 1
            if rate < s['best'][0]:
                s['second'] = s['best']
                s['best'] = (rate, ref)
            elif s['second'] is None or rate < s['second'][0]:
                s['second'] = (rate, ref)
            if rate < 0.02:
                s['near'].append((rate, ref))

rows = []
for query, s in stats.items():
    best_rate, best_ref = s['best']
    second_rate, second_ref = s['second'] if s['second'] else (float('nan'), 'N/A')
    avg_others = (s['sum'] - best_rate) / (s['count'] - 1) if s['count'] > 1 else float('nan')
    ratio = best_rate / avg_others if avg_others > 0 else float('nan')
    note_map = {ref: r for r, ref in s['near'] if ref != best_ref}
    if (second_ref != 'N/A' and second_ref not in note_map
            and not math.isnan(avg_others) and avg_others > 0
            and second_rate / avg_others < 0.1):
        note_map[second_ref] = second_rate
    near_others = sorted(note_map.items(), key=lambda x: x[1])
    notes = ('{' + ','.join(f'{ref}:{r:.4f}' for ref, r in near_others) + '}') if near_others else ''
    rows.append((best_rate, query, best_ref, second_ref, second_rate, avg_others, ratio, notes))
rows.sort()

with open(summary_file, 'w') as out:
    out.write('\t'.join(['QUERY','BEST_MATCH','BEST_RATE','2ND_MATCH','2ND_RATE','AVG_OTHERS','RATIO','NOTES']) + '\n')
    for best_rate, query, best_ref, second_ref, second_rate, avg_others, ratio, notes in rows:
        out.write(f'{query}\t{best_ref}\t{best_rate:.6f}\t{second_ref}\t{second_rate:.6f}\t{avg_others:.6f}\t{ratio:.6f}\t{notes}\n')

print(f"Summary: {summary_file}  ({len(rows)} queries, "
      f"{sum(1 for r in rows if not math.isnan(r[6]) and r[6] < 0.1)} likely duplicates)")

best_rates   = [r[0] for r in rows]
second_rates = [r[4] for r in rows if not math.isnan(r[4])]
avg_rates    = [r[5] for r in rows if not math.isnan(r[5])]
# r[6]=ratio, r[7]=notes

def kde_xy(data, n=500):
    if len(data) < 2:
        return None, None
    kde = gaussian_kde(data)
    xs = np.linspace(min(data), max(data), n)
    return xs, kde(xs)

fig, ax = plt.subplots(figsize=(10, 5))
layers = [
    (avg_rates,    'Avg others',     '#b0b0b0', 1.2, 0.40, '--', 1),
    (second_rates, '2nd best match', '#f5a42a', 2.0, 0.65, '-',  2),
    (best_rates,   'Best match',     '#1a6db5', 2.8, 1.00, '-',  3),
]
for data, label, color, lw, alpha, ls, zorder in layers:
    xs, ys = kde_xy(data)
    if xs is None:
        continue
    ax.plot(xs, ys, color=color, linewidth=lw, alpha=alpha,
            label=label, linestyle=ls, zorder=zorder)
    ax.fill_between(xs, ys, alpha=alpha * 0.25, color=color, zorder=zorder)

ax.set_xlabel('Discordance rate')
ax.set_ylabel('Density')
ax.set_title(f'Gtcheck discordance rate distributions\n{prefix}', fontsize=11)
ax.legend()
plt.tight_layout()
plt.savefig(out_png, dpi=150)
print(f"Plot:    {out_png}")
PYEOF
  >>>

  output {
    File gtcheck_summary = prefix + ".summary.tsv"
    File gtcheck_plot    = prefix + ".distribution.png"
  }
}
