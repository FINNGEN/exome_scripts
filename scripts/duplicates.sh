#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: $0 <plink_prefix1> <plink_prefix2> <output_prefix>"
    echo ""
    echo "  plink_prefix1 : First plink dataset (without .bed extension)"
    echo "  plink_prefix2 : Second plink dataset (without .bed extension)"
    echo "  output_prefix : Prefix for output files"
    echo ""
    echo "Runs KING --duplicate on two datasets to identify cross-dataset"
    echo "sample pairs with kinship > 0.354 (duplicates / MZ twins)."
    echo ""
    echo "Output: <output_prefix>.con  — cross-duplicate pairs"
    exit 1
}

[[ $# -ne 3 ]] && usage

DATASET1="$1"
DATASET2="$2"
OUTPUT="$3"

[[ ! -f "${DATASET1}.bed" ]] && { echo "ERROR: ${DATASET1}.bed not found"; exit 1; }
[[ ! -f "${DATASET2}.bed" ]] && { echo "ERROR: ${DATASET2}.bed not found"; exit 1; }

echo "Dataset 1: $DATASET1 ($(wc -l < "${DATASET1}.fam") samples, $(wc -l < "${DATASET1}.bim") SNPs)"
echo "Dataset 2: $DATASET2 ($(wc -l < "${DATASET2}.fam") samples, $(wc -l < "${DATASET2}.bim") SNPs)"
echo ""

# KING with two comma-separated bed files computes kinship across both datasets
echo "Running KING --duplicate..."
king -b "${DATASET1}.bed","${DATASET2}.bed" --duplicate --prefix "$OUTPUT"

# Filter to cross-dataset pairs only: keep pairs where IID1 and IID2 are from different datasets
if [[ -f "${OUTPUT}.con" ]]; then
    awk 'NR==FNR { ids[$1]=1; next }
         FNR==1  { print; next }
         { if (($2 in ids) != ($4 in ids)) print }
    ' <(awk '{print $2}' "${DATASET1}.fam") "${OUTPUT}.con" > "${OUTPUT}.con.tmp"
    mv "${OUTPUT}.con.tmp" "${OUTPUT}.con"
fi

echo ""
echo "Done."
if [[ -f "${OUTPUT}.con" ]]; then
    N_DUPS=$(tail -n +2 "${OUTPUT}.con" | wc -l)
    echo "  Duplicate pairs: ${N_DUPS}"
    echo "  Output:          ${OUTPUT}.con"
    if [[ $N_DUPS -gt 0 ]]; then
        echo ""
        echo "First few duplicate pairs:"
        head -5 "${OUTPUT}.con"
    fi
else
    echo "  No .con file produced — no duplicates found or KING reported no results."
fi
