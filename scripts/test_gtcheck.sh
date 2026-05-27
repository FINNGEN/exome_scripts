#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: $0 <query_vcf> <ref_vcf> [output_prefix] [--parallel N]"
    echo ""
    echo "  query_vcf    : Query VCF/BCF file (must be indexed)"
    echo "  ref_vcf      : Reference VCF/BCF file (must be indexed)"
    echo "  output_prefix: Output prefix (default: gtcheck_out)"
    echo "  --parallel N : Chunk size for parallel runs (default: 1000000)"
    exit 1
}

QUERY_VCF=""
REF_VCF=""
PREFIX="gtcheck_out"
CHUNK_SIZE=1000000

while [[ $# -gt 0 ]]; do
    case "$1" in
        --parallel) CHUNK_SIZE="$2"; shift 2 ;;
        -*)         echo "ERROR: Unknown flag $1"; usage ;;
        *)
            if   [[ -z "$QUERY_VCF" ]]; then QUERY_VCF="$1"
            elif [[ -z "$REF_VCF"   ]]; then REF_VCF="$1"
            else                              PREFIX="$1"
            fi
            shift ;;
    esac
done

[[ -z "$QUERY_VCF" || -z "$REF_VCF" ]] && usage

[[ ! -f "$QUERY_VCF" ]] && { echo "ERROR: VCF not found: $QUERY_VCF"; exit 1; }
[[ ! -f "$REF_VCF"   ]] && { echo "ERROR: VCF not found: $REF_VCF";   exit 1; }

OUTPUT="${PREFIX}.gtcheck"
SUMMARY="${PREFIX}.summary.tsv"

echo "Query: $QUERY_VCF"
echo "Ref:   $REF_VCF"
echo "Chunks: $CHUNK_SIZE samples each"
echo ""

# --- Unphase + split + gtcheck -----------------------------------------------
if [[ -f "$OUTPUT" ]]; then
    echo "Skipping gtcheck — $OUTPUT already exists"
else
    OUT_DIR=$(dirname "$PREFIX")
    REF_NAME=$(basename "$REF_VCF" | sed 's/\.vcf\.gz$//;s/\.bcf$//')
    REF_VCF_UNPHASED="${OUT_DIR}/${REF_NAME}_unphased.vcf.gz"
    QUERY_UNPHASED="${PREFIX}_query_unphased.vcf.gz"

    if [[ ! -f "$REF_VCF_UNPHASED" ]] || [[ ! -f "${REF_VCF_UNPHASED}.tbi" ]]; then
        echo "Unphasing ref VCF..."
        bcftools view "$REF_VCF" | sed 's/\([0-9]\)|\([0-9]\)/\1\/\2/g' | bcftools view -O z --threads 4 -o "$REF_VCF_UNPHASED"
        bcftools index --tbi --threads 4 "$REF_VCF_UNPHASED"
    else
        echo "Reusing existing $REF_VCF_UNPHASED"
    fi

    if [[ ! -f "$QUERY_UNPHASED" ]] || [[ ! -f "${QUERY_UNPHASED}.tbi" ]]; then
        echo "Unphasing query VCF..."
        bcftools view "$QUERY_VCF" | sed 's/\([0-9]\)|\([0-9]\)/\1\/\2/g' | bcftools view -O z --threads 4 -o "$QUERY_UNPHASED"
        bcftools index --tbi --threads 4 "$QUERY_UNPHASED"
    else
        echo "Reusing existing $QUERY_UNPHASED"
    fi

    # Split sample IDs into chunks for parallel subsetting
    bcftools query -l "$REF_VCF_UNPHASED" \
        | split -d -l "$CHUNK_SIZE" - "${PREFIX}_chunk_"
    mapfile -t CHUNKS < <(ls "${PREFIX}_chunk_"*)

    # One pipeline per chunk: subset+index → gtcheck → cleanup
    PIPELINE_SH="${PREFIX}_pipeline.sh"
    for chunk in "${CHUNKS[@]}"; do
        echo "bcftools view -S ${chunk} -O b --write-index -o ${chunk}_ref.bcf ${REF_VCF_UNPHASED} && /usr/bin/time -v bcftools gtcheck --no-HWE-prob -g ${chunk}_ref.bcf ${QUERY_UNPHASED} > ${chunk}.gtcheck 2> ${chunk}.memlog && rm ${chunk}_ref.bcf ${chunk}_ref.bcf.csi"
    done > "$PIPELINE_SH"
    parallel -j "$(nproc)" < "$PIPELINE_SH"

    printf "\n%-25s %10s %12s\n" "CHUNK" "SAMPLES" "PEAK_RSS_MB"
    for chunk in "${CHUNKS[@]}"; do
        n=$(wc -l < "$chunk")
        peak_kb=$(grep "Maximum resident set size" "${chunk}.memlog" 2>/dev/null | awk '{print $NF}')
        printf "%-25s %10d %12d\n" "$(basename "$chunk")" "$n" "$(( ${peak_kb:-0} / 1024 ))"
    done
    echo ""

    cat "${PREFIX}_chunk_"*.gtcheck > "$OUTPUT"
    rm -f "${PREFIX}_chunk_"* "$PIPELINE_SH"
fi
echo "Output: $OUTPUT"
echo ""

# --- Summary --------------------------------------------------------------
python3 - "$OUTPUT" << 'EOF' > "$SUMMARY"
import sys, math

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

header = (f"{'QUERY':<30} {'BEST_MATCH':<30} {'BEST_RATE':>10}"
          f" {'2ND_MATCH':<30} {'2ND_RATE':>10} {'AVG_OTHERS':>12} {'RATIO':>8}  NOTES")
print(header)

rows = []
for query, s in stats.items():
    best_rate, best_ref     = s['best']
    second_rate, second_ref = s['second'] if s['second'] else (float('nan'), 'N/A')
    avg_others = (s['sum'] - best_rate) / (s['count'] - 1) if s['count'] > 1 else float('nan')
    ratio = best_rate / avg_others if avg_others > 0 else float('nan')
    # near_others: absolute rate < 2% (excluding best) + 2nd match if its ratio is also low
    note_map = {ref: r for r, ref in s['near'] if ref != best_ref}
    if (second_ref != 'N/A' and second_ref not in note_map
            and not math.isnan(avg_others) and avg_others > 0
            and second_rate / avg_others < 0.1):
        note_map[second_ref] = second_rate
    near_others = sorted(note_map.items(), key=lambda x: x[1])
    notes = ('{' + ','.join(f'{ref}:{r:.4f}' for ref, r in near_others) + '}') if near_others else ''
    rows.append((best_rate, query, best_ref, second_ref, second_rate, avg_others, ratio, notes))

rows.sort(key=lambda x: x[0])
for best_rate, query, best_ref, second_ref, second_rate, avg_others, ratio, notes in rows:
    print(f"{query:<30} {best_ref:<30} {best_rate:>10.4f}"
          f" {second_ref:<30} {second_rate:>10.4f} {avg_others:>12.4f} {ratio:>8.4f}  {notes}")
EOF

# --- Plot distributions ---------------------------------------------------
python3 - "$SUMMARY" "${PREFIX}.distribution.png" "$PREFIX" << 'EOF'
import sys
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

summary_file, out_png, prefix = sys.argv[1], sys.argv[2], sys.argv[3]

best_rates, second_rates, avg_rates = [], [], []
with open(summary_file) as f:
    for i, line in enumerate(f):
        if i == 0:
            continue
        parts = line.split()
        if len(parts) < 6:
            continue
        try:
            best_rates.append(float(parts[2]))
            second_rates.append(float(parts[4]))
            avg_rates.append(float(parts[5]))
        except (ValueError, IndexError):
            continue

if not best_rates:
    print("No data to plot", file=sys.stderr)
    sys.exit(0)

def kde_xy(data, n=500):
    clean = [x for x in data if not math.isnan(x)]
    if len(clean) < 2:
        return None, None
    kde = gaussian_kde(clean)
    xs = np.linspace(min(clean), max(clean), n)
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
print(f"Plot saved: {out_png}")
EOF
