#!/bin/bash


if [[ $# -lt 2 ]]; then
  echo "Usage: $0 input.vcf.gz N"
  exit 1
fi

input_vcf="$1"
N="$2"

test_vcf="benchmark_subset.vcf.gz"
regions_dir="benchmark_regions"
chunks=$(nproc)
out_native="benchmark_result_native.vcf.gz"
out_split_dir="benchmark_split"
out_split="benchmark_result_split.vcf.gz"

# 1. Prepare subset of N variants
echo "Building test VCF with first $N variants + header..."
bcftools view -h "$input_vcf" > benchmark_subset.vcf
bcftools view -H "$input_vcf" | head -n "$N" >> benchmark_subset.vcf
bgzip -f benchmark_subset.vcf
tabix -p vcf "$test_vcf"

# 2. Native bcftools view with threads
echo "Timing native bcftools view --threads=$chunks ..."
/usr/bin/time -f "Native method: %E real, %U user, %S sys" \
  bash -c "bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' '$test_vcf' -Ou \
  | bcftools view -e 'FILTER~\"NO_HQ_GENOTYPES\"' -Oz --threads $chunks -o '$out_native'" || true
tabix -p vcf "$out_native"

# 3. Split by regions and parallel
rm -rf "$regions_dir" "$out_split_dir"
mkdir -p "$regions_dir" "$out_split_dir"
echo "Extracting split regions..."
bcftools query -f'%CHROM\t%POS\t%POS\n' "$test_vcf" > "$regions_dir/all_positions.txt"
split -n l/$chunks -d "$regions_dir/all_positions.txt" "$regions_dir/chunk_"
for f in "$regions_dir"/chunk_*; do
  if [[ -s "$f" ]]; then
    awk 'NR==1{chrom=$1;start=$2} {end=$2} END{if(NR>0) print chrom"\t"start"\t"end}' "$f" > "$regions_dir/region_$(basename "$f")"
  fi
done

echo "Timing region split + GNU parallel ..."
# shellcheck disable=SC2016
/usr/bin/time -f "Split/parallel method: %E real, %U user, %S sys" \
  bash -c "
    ls $regions_dir/region_chunk_* | parallel -j $chunks 'bcftools view -R {} \"$test_vcf\" -Ou |
      bcftools view -e \"FILTER~\\\"NO_HQ_GENOTYPES\\\"\" -Ou | 
      bcftools annotate --set-id +\"%CHROM\_%POS\_%REF\_%ALT\" -Oz -o $out_split_dir/chunk_{/.}.vcf.gz'
  " || true
for chunk in "$out_split_dir"/chunk_*.vcf.gz; do
  [[ -f "$chunk" ]] && tabix -p vcf "$chunk"
done

# 4. Concatenate chunks
echo "Concatenating region-split chunks to $out_split ..."
if ls "$out_split_dir"/chunk_*.vcf.gz 1> /dev/null 2>&1; then
  bcftools concat -n -Oz -o "$out_split" "$out_split_dir"/chunk_*.vcf.gz
  tabix -p vcf "$out_split"
else
  echo "Warning: No chunks found to concatenate"
fi

echo ""
echo "=== Results ==="
[[ -f "$out_native" ]] && echo "Native bcftools output: $out_native ($(ls -lh "$out_native" | awk '{print $5}'))"
[[ -f "$out_split" ]] && echo "Region-split output: $out_split ($(ls -lh "$out_split" | awk '{print $5}'))"
