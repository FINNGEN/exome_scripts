#!/usr/bin/env bash
set -euo

usage() {
    echo "Usage: $0 <vcf_or_plink_input> <prefix> <snp_list> [--afreq <afreq_file>]"
    echo ""
    echo "  vcf_or_plink_input : Input VCF/BCF file, or plink prefix (if .bed exists)"
    echo "  prefix             : Output file prefix; also prepended to sample IIDs (PREFIX_OLDIID)"
    echo "  snp_list           : File with SNP IDs to extract (one per line)"
    echo "  --afreq            : Optional .afreq file to speed up plink2 (plink input only)"
    echo ""
    echo "Output: <prefix>.bed/.bim/.fam with --maj-ref force and IIDs renamed to PREFIX_OLDIID"
    exit 1
}

[[ $# -lt 3 ]] && usage

INPUT="$1"
PREFIX="$2"
SNP_LIST="$3"
AFREQ="${INPUT%.bed}.afreq"
shift 3

while [[ $# -gt 0 ]]; do
    case "$1" in
        --afreq) AFREQ="$2"; shift 2 ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# normalize: if user passed a plink prefix (without .bed), resolve to .bed path
PLINK_PREFIX="${INPUT%.bed}"
if [[ -f "${PLINK_PREFIX}.bed" ]]; then
    INPUT="${PLINK_PREFIX}.bed"
    echo "Input:    plink bfile ${PLINK_PREFIX}"
else
    echo "Input:    VCF $INPUT"
    [[ ! -f "$INPUT" ]] && { echo "ERROR: VCF not found: $INPUT"; exit 1; }
fi

[[ ! -f "$SNP_LIST" ]] && { echo "ERROR: SNP list not found: $SNP_LIST"; exit 1; }

THREADS=$(nproc)

echo "Prefix:   $PREFIX"
[[ -n "$AFREQ" ]] && echo "Afreq:    $AFREQ"
echo "SNP list: $SNP_LIST ($(wc -l < "$SNP_LIST") SNPs)"
echo ""

if [[ "$INPUT" == *.bed ]]; then
  echo "=== Plink to Plink: $PREFIX ==="
  INPUT_FLAGS="--bfile ${INPUT%.bed}${AFREQ:+ --read-freq $AFREQ}"
else
  echo "=== VCF to Plink: $PREFIX ==="
  INPUT_FLAGS="--vcf $INPUT --double-id --max-alleles 2"
fi

plink2 \
  $INPUT_FLAGS \
  --extract "$SNP_LIST" \
  --autosome \
  --maj-ref force \
  --make-bed \
  --out "$PREFIX" \
  --threads $THREADS \
  --allow-extra-chr

echo ""
echo "Renaming IIDs to ${PREFIX}_OLDIID..."
awk -v prefix="$PREFIX" 'BEGIN{OFS="\t"} {print $1, $2, $1, prefix "_" $2}' "${PREFIX}.fam" > id_mapping.txt
plink2 \
  --bfile "$PREFIX" \
  --update-ids id_mapping.txt \
  --make-just-fam \
  --out "$PREFIX"
rm -f id_mapping.txt

echo ""
echo "Done."
echo "  SNPs:    $(wc -l < "${PREFIX}.bim")"
echo "  Samples: $(wc -l < "${PREFIX}.fam")"
