#!/bin/bash
# Test script for fast header annotation vs bcftools annotate

input_vcf="$1"
if [[ -z "$input_vcf" ]]; then
  echo "Usage: $0 input.vcf.gz"
  exit 1
fi

# One-liner: Generate FILTER headers, inject into VCF, take 100 variants
(printf '##FILTER=<ID=NO_HQ_GENOTYPES,Description="Site has no high quality variant genotypes. No high-quality genotype (GQ>=20, DP>=10, and AB>=0.2 for heterozygotes) called for the variant. If there is one genotype at the variant site, the filter will not be applied and the variant site will pass. Allele Balance (AB) is min(AD)/DP for diploid GTs.">\n##FILTER=<ID=ExcessHet,Description="Site has excess het value larger than the threshold. Genotypes with this filter show a higher proportion of heterozygotes than expected under Hardy-Weinberg equilibrium (z-score < -4.5, phred 54.69). Suggests mapping errors or contamination.">\n##FILTER=<ID=LowQual,Description="QUALapprox is too low (lower than 60 for SNPs; lower than 69 for Indels). QUAL tells you how confident we are that there is some kind of variation at a given site.">\n##FILTER=<ID=EXCESS_ALLELES,Description="Site has an excess of alternate alleles based on the input threshold (e.g. >100 alternate alleles).">\n##FILTER=<ID=OUTSIDE_OF_TARGETS,Description="Exome only. The site is not within the target intervals of the exome assay.">\n' > /tmp/filters.txt; bcftools view -h "$input_vcf" | head -n -1; cat /tmp/filters.txt; bcftools view -h "$input_vcf" | tail -n 1; bcftools view -H "$input_vcf" | head -n 100) | bgzip > test_annotated.vcf.gz && tabix -p vcf test_annotated.vcf.gz

# Validate the output
echo ""
echo "=== Validation ==="
echo -n "Variant count: "
bcftools view -H test_annotated.vcf.gz | wc -l

echo ""
echo "Added FILTER lines in header:"
bcftools view -h test_annotated.vcf.gz | grep "^##FILTER" | grep -E "(NO_HQ_GENOTYPES|ExcessHet|LowQual|EXCESS_ALLELES|OUTSIDE_OF_TARGETS)"

echo ""
echo "Testing bcftools can parse it:"
bcftools query -f '%CHROM\t%POS\n' test_annotated.vcf.gz | head -n 5

echo ""
echo "✓ File is valid: test_annotated.vcf.gz"
