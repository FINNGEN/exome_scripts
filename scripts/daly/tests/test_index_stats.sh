#!/bin/bash
# Validate that bcftools index -s gives us accurate last position

input_vcf="$1"
if [[ -z "$input_vcf" ]]; then
  echo "Usage: $0 input.vcf.gz"
  exit 1
fi

echo "=== Testing index stats vs actual variant positions ==="
echo ""

# Method 1: Current slow method (zcat entire file)
echo "Method 1: zcat and awk entire VCF..."
time {
  read chrom1 first1 last1 < <(zcat "$input_vcf" | awk '!/^#/ {chrom=$1; if(!first) first=$2; last=$2} END{print chrom, first, last}')
}
echo "  Chrom: $chrom1"
echo "  First: $first1"
echo "  Last:  $last1"
echo ""

# Method 2: Fast method using index
echo "Method 2: Using bcftools index -s..."
time {
  read chrom2 last2 < <(bcftools index -s "$input_vcf" | awk '{print $1, $2}')
  first2=$(bcftools query -f'%POS\n' "$input_vcf" | head -n 1)
}
echo "  Chrom: $chrom2"
echo "  First: $first2"
echo "  Last:  $last2"
echo ""

# Method 3: Just verify - get actual last variant position
echo "Method 3: Getting actual last variant position (zcat + tail)..."
actual_last=$(zcat "$input_vcf" | awk '!/^#/ {last=$2} END{print last}')
echo "  Actual last variant position: $actual_last"
echo ""

# Comparison
echo "=== Comparison ==="
echo "First position match: $([ "$first1" == "$first2" ] && echo "✓ YES" || echo "✗ NO ($first1 vs $first2)")"
echo "Last position (scan vs index): $([ "$last1" == "$last2" ] && echo "✓ MATCH" || echo "✗ DIFFER ($last1 vs $last2)")"
echo "Last position (index vs actual): $([ "$last2" == "$actual_last" ] && echo "✓ MATCH" || echo "✗ DIFFER ($last2 vs $actual_last)")"
echo ""

if [[ "$last2" == "$actual_last" ]]; then
  echo "✓ SUCCESS: bcftools index -s column 2 gives exact last variant position"
  echo "  You can safely use the fast method!"
else
  echo "✗ WARNING: bcftools index -s column 2 ($last2) != actual last variant ($actual_last)"
  echo "  Difference: $((last2 - actual_last)) bp"
fi
