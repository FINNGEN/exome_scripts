#!/bin/bash
# Filter VCF: set low-quality genotypes to missing and recalculate AC
# Usage: ./filter_vcf.sh [input.vcf.gz] [threads]

input_file="${1:-/mnt/disks/data/exome/test/test/ADPKD_test.vcf.gz}"
threads="${2:-$(nproc)}"
output_file="/mnt/disks/data/exome/test/test/filtered.vcf.gz"

echo "Input:   $input_file"
echo "Output:  $output_file"
echo "Threads: $threads"
echo ""

echo "Executing command:"
echo "---"
cat << 'EOF'
bcftools +setGT "$input_file" -Ou -- -t q -n . -i "FORMAT/DP<10 | FORMAT/GQ<20" | \
  bcftools +fill-tags -Ou -- -t AC | \
  bcftools view -i 'AC>0 & ALT!="*"' --threads "$threads" -Ou | \
  bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' --threads "$threads" -Oz -o "$output_file"
EOF
echo "---"
echo ""

bcftools +setGT "$input_file" -Ou -- -t q -n . -i "FORMAT/DP<10 | FORMAT/GQ<20" | \
  bcftools +fill-tags -Ou -- -t AC | \
  bcftools view -i 'AC>0 & ALT!="*"' --threads "$threads" -Ou | \
  bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' --threads "$threads" -Oz -o "$output_file"

echo ""
echo "Indexing..."
tabix -p vcf "$output_file"

echo ""
echo "Done! Output: $output_file"
