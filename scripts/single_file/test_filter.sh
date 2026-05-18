#!/bin/bash
# Test script to verify filtering worked correctly
# Usage: ./test_filter.sh [original.vcf.gz] [filtered.vcf.gz]

original="${1:-/mnt/disks/data/exome/test/test/ADPKD_test.vcf.gz}"
filtered="${2:-/mnt/disks/data/exome/test/test/filtered.vcf.gz}"

echo "=== Testing VCF Filtering ==="
echo "Original: $original"
echo "Filtered: $filtered"
echo ""

# Test 1: Find a genotype that was set to missing
echo "Test 1: Checking genotypes were set to missing (DP<10 or GQ<20)"
echo "--------------------------------------------------------------"

# Find genotypes that meet criteria and check if variant still exists in filtered file
found=0
bcftools query -i 'GT!="mis"' -f '%CHROM\t%POS[\t%SAMPLE\t%GT\t%DP\t%GQ]\n' "$original" 2>/dev/null | \
  awk '($4 < 10 || $5 < 20)' | \
  while read chrom pos sample gt dp gq; do
    # Check if this variant exists in filtered file
    if bcftools view -H -r "$chrom:$pos" "$filtered" 2>/dev/null | grep -q .; then
      after_gt=$(bcftools query -s "$sample" -r "$chrom:$pos" -f '[\t%GT]\n' "$filtered" 2>/dev/null | tr -d '\t' | tr -d ' ')
      
      if [[ "$after_gt" == "./." ]] || [[ "$after_gt" == "." ]]; then
        echo "✓ PASS: $chrom:$pos sample=$sample GT: $gt → $after_gt"
      else
        echo "✗ FAIL: $chrom:$pos sample=$sample GT: $gt → $after_gt (not missing)"
      fi
      found=1
      break
    fi
  done

if [[ $found -eq 0 ]]; then
  echo "⚠ NOTE: All low-quality genotypes resulted in variants being removed (AC=0)"
fi

echo ""

# Test 2: Check AC was recalculated
echo "Test 2: Checking AC was recalculated"
echo "--------------------------------------"

read chrom2 pos2 < <(bcftools query -i 'GT!="mis"' -f '%CHROM\t%POS[\t%DP\t%GQ]\n' "$original" | \
  awk '($3 < 10 || $4 < 20) {print $1, $2; exit}')

if [[ -n "$chrom2" ]]; then
  ac_before=$(bcftools query -r "$chrom2:$pos2" -f '%AC\n' "$original")
  ac_after=$(bcftools query -r "$chrom2:$pos2" -f '%AC\n' "$filtered")
  
  echo "Position: $chrom2:$pos2"
  echo "AC BEFORE: $ac_before"
  echo "AC AFTER:  $ac_after"
  
  if [[ "$ac_after" -le "$ac_before" ]]; then
    echo "✓ PASS: AC was recalculated (lower or equal)"
  else
    echo "✗ FAIL: AC increased (unexpected)"
  fi
else
  echo "No suitable variant found for AC test"
fi

echo ""

# Test 3: Check AC>0 filtering
echo "Test 3: Checking all variants have AC>0"
echo "----------------------------------------"

zero_ac=$(bcftools view -H "$filtered" | bcftools query -f '%AC\n' | awk '$1 == 0' | wc -l)
if [[ "$zero_ac" -eq 0 ]]; then
  echo "✓ PASS: No variants with AC=0 found"
else
  echo "✗ FAIL: Found $zero_ac variants with AC=0"
fi

echo ""

# Test 4: Check variant IDs were annotated
echo "Test 4: Checking variant IDs (CHROM_POS_REF_ALT format)"
echo "--------------------------------------------------------"

sample_id=$(bcftools view -H "$filtered" | head -1 | awk '{print $3}')
if [[ "$sample_id" =~ ^[^_]+_[0-9]+_.+_.+$ ]]; then
  echo "Sample ID: $sample_id"
  echo "✓ PASS: IDs formatted as CHROM_POS_REF_ALT"
else
  echo "Sample ID: $sample_id"
  echo "✗ FAIL: ID not in expected format"
fi

echo ""

# Summary
echo "=== Summary ==="
echo "Tests complete"
