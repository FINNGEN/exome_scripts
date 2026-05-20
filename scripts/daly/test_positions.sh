#!/bin/bash
# Usage: test_positions.sh <vcf.gz>
# Tests to find the actual first and last variant positions in a VCF
# Requires: bcftools

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <vcf.gz>"
  exit 1
fi

VCF="$1"

echo -e "contig\tlength\tvariants\tfirst_pos\tlast_pos"

# For each contig using bcftools index -s
bcftools index -s "$VCF" | while read -r chrom length nvar; do
  if [[ "$nvar" == "0" ]]; then
    echo -e "$chrom\t$length\t0\tNA\tNA"
    continue
  fi

  # Get first position from bcftools view -H | head -n1
  first_pos=$(bcftools view -H -r "$chrom" "$VCF" | head -n 1 | awk '{print $2}')
  echo "[$chrom] First position: $first_pos" >&2

  # Binary search to find the last position with variants
  low=$first_pos
  high=$length
  last_pos=""
  
  echo "[$chrom] Starting binary search from $low to $high" >&2
  
  while (( low <= high )); do
    mid=$(( (low + high) / 2 ))
    
    echo "[$chrom] Testing range $low-$high (midpoint: $mid)" >&2
    
    # Quick check: is there ANY variant in the upper half? (exit as soon as we find one)
    has_variant=$(bcftools view -H -r "$chrom:$mid-$high" "$VCF" 2>/dev/null | head -n 1)
    
    if [[ -n "$has_variant" ]]; then
      # Found at least one variant in upper half, search higher
      echo "[$chrom] Found variant(s) in range $mid-$high, searching higher" >&2
      low=$(( mid + 1 ))
    else
      # No variants in upper half, search lower
      echo "[$chrom] No variants in range $mid-$high, searching lower" >&2
      high=$(( mid - 1 ))
    fi
  done
  
  # At this point, 'high' is the highest position where variants exist
  # Do a very narrow search to get the exact last variant (just last 10kb)
  if (( high >= first_pos )); then
    echo "[$chrom] Binary search narrowed to position $high, finding exact last variant" >&2
    # Search only the last 10kb to quickly find the last variant
    search_start=$(( high > 10000 ? high - 10000 : first_pos ))
    last_pos=$(bcftools view -H -r "$chrom:$search_start-$high" "$VCF" | tail -n 1 | awk '{print $2}')
    echo "[$chrom] Last variant at position $last_pos" >&2
  fi

  # Fallback if still not found
  if [[ -z "$last_pos" ]]; then
    last_pos=$(bcftools view -H -r "$chrom" "$VCF" | tail -n 1 | awk '{print $2}')
    [[ -z "$last_pos" ]] && last_pos="NA"
  fi

  echo -e "$chrom\t$length\t$nvar\t$first_pos\t$last_pos"
done
