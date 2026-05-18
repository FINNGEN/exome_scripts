#!/bin/bash
# Parallel VCF filtering by chromosome
# Sets low-quality genotypes to missing, recalculates AC, filters variants, and annotates IDs
# Usage: ./parallel_filter_by_chrom.sh <input.vcf.gz> [--chroms 1 2 3 ...] [--cpus N]

usage() {
  echo "Usage: $0 <input.vcf.gz> [--chroms chr1 chr2 ...] [--cpus N]"
  echo ""
  echo "This script:"
  echo "  1. Sets genotypes to missing where FORMAT/DP<10 | FORMAT/GQ<20"
  echo "  2. Recalculates AC (allele count)"
  echo "  3. Filters to keep only variants with AC>0 & ALT!=\"*\""
  echo "  4. Annotates variant IDs as CHROM_POS_REF_ALT"
  echo ""
  echo "  --chroms  Space-separated list of chromosomes to process (default: all)"
  echo "  --cpus    Number of parallel jobs (default: nproc)"
  echo ""
  echo "Processing is parallelized by chromosome."
  exit 1
}

[[ $# -lt 1 ]] && usage

input_file="$1"
shift

CHUNKS=$(nproc)
NORM_FASTA="/mnt/disks/data/exome/fasta/Homo_sapiens_assembly38.fasta"
CHROMS_ARG=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --chroms)
      shift
      while [[ $# -gt 0 && "$1" != --* ]]; do
        CHROMS_ARG+=("$1")
        shift
      done
      ;;
    --cpus)
      CHUNKS="$2"
      shift 2
      ;;
    *)
      echo "Unknown argument: $1"
      usage
      ;;
  esac
done

BASENAME=$(basename "$input_file")
BASENAME="${BASENAME%.vcf.gz}"
BASENAME="${BASENAME%.vcf.bgz}"
BASENAME="${BASENAME%.bcf.gz}"
BASENAME="${BASENAME%.bcf}"
OUTPUT_FILE="${BASENAME}_QC_ANNOTATED.vcf.gz"

echo "=== Parallel Filter by Chromosome: Starting ==="
echo "Input file:  $input_file"
echo "Output file: $OUTPUT_FILE"
echo "CPU cores:   $CHUNKS"

# Check if input file exists
if [[ ! -f "$input_file" ]]; then
  echo "Error: Input file not found: $input_file"
  exit 1
fi

# Check if index exists, create if not
if [[ ! -f "$input_file.tbi" ]] && [[ ! -f "$input_file.csi" ]]; then
  echo "Creating index..."
  tabix -p vcf "$input_file"
fi

# Get list of chromosomes — use --chroms if provided, otherwise read from VCF index
if [[ ${#CHROMS_ARG[@]} -gt 0 ]]; then
  chromosomes=("${CHROMS_ARG[@]}")
  echo "Using specified chromosomes: ${chromosomes[*]}"
else
  echo "Extracting chromosome list from VCF index..."
  mapfile -t chromosomes < <(bcftools index -s "$input_file" | awk '{print $1}')
  echo "Found ${#chromosomes[@]} chromosomes: ${chromosomes[*]}"
fi

num_chroms=${#chromosomes[@]}

if [[ $num_chroms -eq 0 ]]; then
  echo "Error: No chromosomes found"
  exit 1
fi

# Create chromosome list file for parallel processing
rm -f chrom_list.txt
for chrom in "${chromosomes[@]}"; do
  echo "$chrom" >> chrom_list.txt
done

echo "Processing $num_chroms chromosomes in parallel (max $CHUNKS jobs)..."

# Detect chromosome naming: if VCF uses non-chr names (e.g. "1" not "chr1"),
# permanently rename to chr prefix so output is always chr-prefixed
FIRST_CHROM=$(bcftools index -s "$input_file" | awk 'NR==1{print $1}')
RENAME_TO_CHR=""
if [[ "$FIRST_CHROM" != chr* ]]; then
  echo "VCF uses non-chr chromosome names — renaming to chr prefix in output"
  RENAME_TO_CHR=$(mktemp)
  for i in $(seq 1 22) X Y MT M; do
    echo "$i chr$i" >> "$RENAME_TO_CHR"
  done
fi

# Generate a per-chromosome script and run them in parallel
# All variables expand at generation time — no function exports needed
SCRIPT_DIR=$(mktemp -d)
for chrom in "${chromosomes[@]}"; do
  output="chunk_${chrom}.vcf.gz"
  if [[ -n "$RENAME_TO_CHR" ]]; then
    cat > "${SCRIPT_DIR}/run_${chrom}.sh" << SCRIPT
#!/bin/bash
echo "Processing chromosome: ${chrom}"
bcftools view -r "${chrom}" "${input_file}" | \\
    tr -d '\0' | \\
    bcftools annotate --rename-chrs "${RENAME_TO_CHR}" -Ou | \\
    bcftools norm -f "${NORM_FASTA}" -m -any -c x -Ou | \\
    bcftools +setGT -Ou -- -t q -n . -i 'FORMAT/DP<10 | FORMAT/GQ<20' | \\
    bcftools +fill-tags -Ou -- -t AC | \\
    bcftools view -i 'AC>0 & ALT!="*"' -Ou | \\
    bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' -Oz -o "${output}"
echo "Completed chromosome: ${chrom}"
SCRIPT
  else
    cat > "${SCRIPT_DIR}/run_${chrom}.sh" << SCRIPT
#!/bin/bash
echo "Processing chromosome: ${chrom}"
bcftools view -r "${chrom}" "${input_file}" | \\
    tr -d '\0' | \\
    bcftools norm -f "${NORM_FASTA}" -m -any -c x -Ou | \\
    bcftools +setGT -Ou -- -t q -n . -i 'FORMAT/DP<10 | FORMAT/GQ<20' | \\
    bcftools +fill-tags -Ou -- -t AC | \\
    bcftools view -i 'AC>0 & ALT!="*"' -Ou | \\
    bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' -Oz -o "${output}"
echo "Completed chromosome: ${chrom}"
SCRIPT
  fi
done

ls "${SCRIPT_DIR}"/run_*.sh | parallel -j "$CHUNKS" 'bash {}'
rm -rf "$SCRIPT_DIR"

echo "Indexing chromosome chunks..."
ls chunk_*.vcf.gz | parallel -j "$CHUNKS" 'tabix -p vcf {}'

echo "Concatenating chromosomes..."
bcftools concat -n -Oz -o "$OUTPUT_FILE" chunk_*.vcf.gz

echo "Indexing final output..."
tabix -p vcf "$OUTPUT_FILE"

echo ""
echo "Output file: $OUTPUT_FILE"

# Cleanup
# echo "Cleaning up temporary files..."
rm -f chunk_*.vcf.gz chunk_*.vcf.gz.tbi chrom_list.txt
[[ -n "$RENAME_TO_CHR" ]] && rm -f "$RENAME_TO_CHR"

echo "=== Parallel Filter by Chromosome: Complete! ==="
