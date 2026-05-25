#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: $0 <query_vcf> <plink_prefix> [output_prefix] [--parallel N]"
    echo ""
    echo "  query_vcf    : Query VCF/BCF file (must be indexed)"
    echo "  plink_prefix : Reference plink dataset (without .bed extension)"
    echo "  output_prefix: Output prefix (default: gtcheck_out)"
    echo "  --parallel N : Chunk size for parallel plink2 runs (default: 1000000)"
    exit 1
}

QUERY_VCF=""
PLINK=""
PREFIX="gtcheck_out"
CHUNK_SIZE=1000000

while [[ $# -gt 0 ]]; do
    case "$1" in
        --parallel) CHUNK_SIZE="$2"; shift 2 ;;
        -*)         echo "ERROR: Unknown flag $1"; usage ;;
        *)
            if   [[ -z "$QUERY_VCF" ]]; then QUERY_VCF="$1"
            elif [[ -z "$PLINK"     ]]; then PLINK="$1"
            else                              PREFIX="$1"
            fi
            shift ;;
    esac
done

[[ -z "$QUERY_VCF" || -z "$PLINK" ]] && usage

[[ ! -f "$QUERY_VCF" ]]   && { echo "ERROR: VCF not found: $QUERY_VCF";  exit 1; }
[[ ! -f "${PLINK}.bed" ]] && { echo "ERROR: ${PLINK}.bed not found";      exit 1; }

PLINK_NAME=$(basename "$PLINK")
RENAMED_FAM="${PREFIX}_${PLINK_NAME}.fam"
OUTPUT="${PREFIX}.gtcheck"
SUMMARY="${PREFIX}.summary.tsv"

echo "Query: $QUERY_VCF"
echo "Ref:   $PLINK"
echo "Chunks: $CHUNK_SIZE samples each"
echo ""

# --- Build renamed fam (skip if already done) -----------------------------
if [[ ! -f "$RENAMED_FAM" ]]; then
    awk -v p="$PLINK_NAME" 'BEGIN{OFS="\t"} {$2 = p"_"$2; print}' "${PLINK}.fam" > "$RENAMED_FAM"
    echo "Renamed fam: $RENAMED_FAM ($(wc -l < "$RENAMED_FAM") samples)"
fi

# --- Plink to VCF + gtcheck -----------------------------------------------
if [[ -f "$OUTPUT" ]]; then
    echo "Skipping gtcheck — $OUTPUT already exists"
else
    awk '{print $1, $2}' "$RENAMED_FAM" | split -d -l "$CHUNK_SIZE" - "${PREFIX}_chunk_"
    mapfile -t CHUNKS < <(ls "${PREFIX}_chunk_"*)

    NJOBS=$(( ${#CHUNKS[@]} < $(nproc) ? ${#CHUNKS[@]} : $(nproc) ))
    TOTAL_MEM_MB=$(awk '/MemAvailable/ {print int($2/1024)}' /proc/meminfo)
    MEM_PER_JOB=$(( TOTAL_MEM_MB / NJOBS ))

    SPLIT_SH="${PREFIX}_split.sh"
    GTCHECK_SH="${PREFIX}_gtcheck.sh"

    # Write split.sh: plink2 + index per chunk
    for chunk in "${CHUNKS[@]}"; do
        echo "plink2 --bfile ${PLINK} --fam ${RENAMED_FAM} --keep ${chunk} --export vcf id-paste=iid bgz --output-chr chrM --out ${chunk}_ref --threads 1 --memory ${MEM_PER_JOB} && bcftools index -t ${chunk}_ref.vcf.gz"
    done > "$SPLIT_SH"
    parallel -j "$(nproc)" < "$SPLIT_SH"

    # Write gtcheck.sh: one gtcheck command per chunk
    for chunk in "${CHUNKS[@]}"; do
        echo "bcftools gtcheck --no-HWE-prob -g ${chunk}_ref.vcf.gz ${QUERY_VCF} > ${chunk}.gtcheck && rm ${chunk}_ref.vcf.gz ${chunk}_ref.vcf.gz.tbi"
    done > "$GTCHECK_SH"
    parallel -j "$(nproc)" < "$GTCHECK_SH"

    cat "${PREFIX}_chunk_"*.gtcheck > "$OUTPUT"
    rm -f "${PREFIX}_chunk_"* "$SPLIT_SH" "$GTCHECK_SH"
fi
echo "Output: $OUTPUT"
echo ""

# --- Summary --------------------------------------------------------------
python3 - "$OUTPUT" << 'EOF' | tee "$SUMMARY"
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

header = (f"{'QUERY':<30} {'BEST_MATCH':<30} {'BEST_RATE':>10}"
          f" {'2ND_MATCH':<30} {'2ND_RATE':>10} {'AVG_OTHERS':>12} {'RATIO':>8}")
print(header)

rows = []
for query, s in stats.items():
    best_rate, best_ref     = s['best']
    second_rate, second_ref = s['second'] if s['second'] else (float('nan'), 'N/A')
    avg_others = (s['sum'] - best_rate) / (s['count'] - 1) if s['count'] > 1 else float('nan')
    ratio = best_rate / avg_others if avg_others > 0 else float('nan')
    rows.append((best_rate, query, best_ref, second_ref, second_rate, avg_others, ratio))

rows.sort(key=lambda x: x[0])
for best_rate, query, best_ref, second_ref, second_rate, avg_others, ratio in rows:
    print(f"{query:<30} {best_ref:<30} {best_rate:>10.4f}"
          f" {second_ref:<30} {second_rate:>10.4f} {avg_others:>12.4f} {ratio:>8.4f}")
EOF
