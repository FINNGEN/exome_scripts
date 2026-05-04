#!/bin/bash
# Parallel VCF filtering by chromosome
# Sets low-quality genotypes to missing, recalculates AC, filters variants, and annotates IDs
# Usage: ./parallel_filter_by_chrom.sh <input.vcf.gz> [cpu_count]

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <input.vcf.gz> [cpu_count]"
  echo ""
  echo "This script:"
  echo "  1. Sets genotypes to missing where FORMAT/DP<10 | FORMAT/GQ<20"
  echo "  2. Recalculates AC (allele count)"
  echo "  3. Filters to keep only variants with AC>0 & ALT!=\"*\""
  echo "  4. Annotates variant IDs as CHROM_POS_REF_ALT"
  echo ""
  echo "Processing is parallelized by chromosome."
  exit 1
fi

input_file="$1"
CHUNKS="${2:-$(nproc)}"

echo "=== Parallel Filter by Chromosome: Starting ==="
echo "Input file: $input_file"
echo "CPU cores: $CHUNKS"

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

# Get list of chromosomes from the VCF index
echo "Extracting chromosome list..."
mapfile -t chromosomes < <(bcftools index -s "$input_file" | awk '{print $1}')
num_chroms=${#chromosomes[@]}

echo "Found $num_chroms chromosomes: ${chromosomes[*]}"

if [[ $num_chroms -eq 0 ]]; then
  echo "Error: No chromosomes found in VCF"
  exit 1
fi

# Create chromosome list file for parallel processing
rm -f chrom_list.txt
for chrom in "${chromosomes[@]}"; do
  echo "$chrom" >> chrom_list.txt
done

echo "Processing $num_chroms chromosomes in parallel (max $CHUNKS jobs)..."

# Process each chromosome in parallel
cat chrom_list.txt | parallel -j "$CHUNKS" '
  chrom={}
  output="chunk_{}.vcf.gz"
  echo "Processing chromosome: $chrom"
  
  bcftools view -r "$chrom" "'"$input_file"'" -Ou | \
    bcftools +setGT -Ou -- -t q -n . -i "FORMAT/DP<10 | FORMAT/GQ<20" | \
    bcftools +fill-tags -Ou -- -t AC | \
    bcftools view -i "AC>0 & ALT!=\"*\"" -Ou | \
    bcftools annotate --set-id +"%CHROM\_%POS\_%REF\_%ALT" -Oz -o "$output"
  
  echo "Completed chromosome: $chrom"
'

echo "Indexing chromosome chunks..."
ls chunk_*.vcf.gz | parallel -j "$CHUNKS" 'tabix -p vcf {}'

echo "Concatenating chromosomes..."
bcftools concat -n -Oz -o filtered.vcf.gz chunk_*.vcf.gz

echo "Indexing final output..."
tabix -p vcf filtered.vcf.gz

echo ""
echo "Output file: filtered.vcf.gz"

# Cleanup
echo "Cleaning up temporary files..."
rm -f chunk_*.vcf.gz chunk_*.vcf.gz.tbi chrom_list.txt

echo "=== Parallel Filter by Chromosome: Complete! ==="
