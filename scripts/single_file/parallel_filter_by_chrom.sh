#!/bin/bash
set -eu

usage() {
  echo "Usage: $0 <input.vcf.gz> [--chroms chr1 chr2 ...] [--cpus N] [--norm-fasta FILE] [--genotype-filter EXPR] [--variant-filter EXPR]"
  echo ""
  echo "  --chroms           Space-separated chromosomes to process (default: all from VCF index)"
  echo "  --cpus             Number of parallel jobs (default: nproc - 1)"
  echo "  --norm-fasta       Reference FASTA (default: GRCh38 on /mnt/disks/data)"
  echo "  --genotype-filter  bcftools filter expression to set GTs missing (default: FORMAT/DP<10 | FORMAT/GQ<20)"
  echo "  --variant-filter   bcftools filter expression to keep variants (default: AC>0 & ALT!=\"*\")"
  exit 1
}

[[ $# -lt 1 || "$1" == --* ]] && usage

INPUT_VCF="$1"
shift

NORM_FASTA="/mnt/disks/data/exome/fasta/Homo_sapiens_assembly38.fasta"
GENOTYPE_FILTER="FORMAT/DP<10 | FORMAT/GQ<20"
VARIANT_FILTER='AC>0 & ALT!="*"'
chromosomes=()
CHUNKS=$(( $(nproc) - 1 ))
if [[ $CHUNKS -lt 1 ]]; then CHUNKS=1; fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    --chroms)
      shift
      while [[ $# -gt 0 && "$1" != --* ]]; do
        chromosomes+=("$1")
        shift
      done
      ;;
    --cpus)             CHUNKS="$2";           shift 2 ;;
    --norm-fasta)       NORM_FASTA="$2";       shift 2 ;;
    --genotype-filter)  GENOTYPE_FILTER="$2";  shift 2 ;;
    --variant-filter)   VARIANT_FILTER="$2";   shift 2 ;;
    *) echo "Unknown argument: $1"; usage ;;
  esac
done

[[ ! -f "$INPUT_VCF" ]] && { echo "ERROR: Input file not found: $INPUT_VCF"; exit 1; }

if [[ ! -f "$INPUT_VCF.tbi" && ! -f "$INPUT_VCF.csi" ]]; then
  echo "Creating index..."
  bcftools index "$INPUT_VCF"
fi

BASENAME=$(basename "$INPUT_VCF")
BASENAME="${BASENAME%.vcf.gz}"
BASENAME="${BASENAME%.vcf.bgz}"
BASENAME="${BASENAME%.bcf.gz}"
BASENAME="${BASENAME%.bcf}"
OUTPUT_VCF="${BASENAME}.QC_ANNOTATED.vcf.gz"

echo "=== Parallel Filter by Chromosome ==="
echo "Input:           $INPUT_VCF"
echo "Output:          $OUTPUT_VCF"
echo "Genotype filter: $GENOTYPE_FILTER"
echo "Variant filter:  $VARIANT_FILTER"
echo "CPU cores:       $CHUNKS"
echo ""

if [[ ${#chromosomes[@]} -eq 0 ]]; then
  mapfile -t chromosomes < <(bcftools index -s "$INPUT_VCF" | awk '{print $1}')
fi
num_chroms=${#chromosomes[@]}
echo "Found $num_chroms chromosomes: ${chromosomes[*]}"

if [[ $num_chroms -eq 0 ]]; then
  echo "Error: No chromosomes found in VCF"
  exit 1
fi

echo "Processing $num_chroms chromosomes in parallel (max $CHUNKS jobs)..."
echo ""

# Detect chromosome naming: if VCF uses non-chr names, rename to chr prefix in output
FIRST_CHROM=$(bcftools index -s "$INPUT_VCF" | awk 'NR==1{print $1}')
RENAME_TO_CHR=""
if [[ "$FIRST_CHROM" != chr* ]]; then
  echo "VCF uses non-chr chromosome names — renaming to chr prefix in output"
  RENAME_TO_CHR=$(mktemp)
  for i in $(seq 1 22) X Y MT M; do
    echo "$i chr$i" >> "$RENAME_TO_CHR"
  done
fi
RENAME_STEP="cat"
[[ -n "$RENAME_TO_CHR" ]] && RENAME_STEP="bcftools annotate --rename-chrs ${RENAME_TO_CHR} -Ou"

# Generate per-chromosome scripts and run in parallel
# Variables expand at generation time — no quoting or function-export issues
SCRIPT_DIR=$(mktemp -d)
for chrom in "${chromosomes[@]}"; do
  safe_chrom=$(echo "$chrom" | sed 's/[*:\/]/_/g')
  output="chunk_${safe_chrom}.vcf.gz"
  cat > "${SCRIPT_DIR}/run_${safe_chrom}.sh" << SCRIPT
#!/bin/bash
echo "Processing chromosome: ${chrom}"
printf '${chrom}\t0\t9999999999\n' > "${output}.region.bed"
bcftools view -R "${output}.region.bed" "${INPUT_VCF}" | \\
  tr -d '\0' | \\
  ${RENAME_STEP} | \\
  bcftools norm -f '${NORM_FASTA}' -m -any -c x -Ou | \\
  bcftools +setGT -Ou -- -t q -n . -i '${GENOTYPE_FILTER}' | \\
  bcftools +fill-tags -Ou -- -t AC | \\
  bcftools view -i '${VARIANT_FILTER}' -Ou | \\
  bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' -Oz -o "${output}"
rm -f "${output}.region.bed"
echo "Completed chromosome: ${chrom}"
SCRIPT
done

ls "${SCRIPT_DIR}"/run_*.sh | parallel -j "$CHUNKS" 'bash {}'
rm -rf "$SCRIPT_DIR"

echo ""
echo "Concatenating chromosomes in original order..."
rm -f chunk_list.txt
for chrom in "${chromosomes[@]}"; do
  safe_chrom=$(echo "$chrom" | sed 's/[*:\/]/_/g')
  echo "chunk_${safe_chrom}.vcf.gz" >> chunk_list.txt
done

bcftools concat -n -Oz -o "$OUTPUT_VCF" -f chunk_list.txt

echo "Indexing final output..."
tabix -p vcf "$OUTPUT_VCF"

# Cleanup
rm -f chunk_*.vcf.gz chunk_list.txt
[[ -n "$RENAME_TO_CHR" ]] && rm -f "$RENAME_TO_CHR"
echo "=== Complete! ==="
