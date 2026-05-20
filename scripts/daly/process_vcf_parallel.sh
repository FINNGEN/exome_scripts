#!/bin/bash

set -euo pipefail

input_file=annotated.fimm-daly_finnish_gvs_bge_callset_1_padded_split_FINBBonly_chr21.vcf.bgz
CHUNKS=$(nproc)

echo "Creating $CHUNKS region files..."
# Get chromosome, first position, and last position
read chrom first_pos last_pos < <(bcftools query -f'%CHROM\t%POS\n' "$input_file" | awk 'NR==1{chrom=$1; first=$2} {last=$2} END{print chrom, first, last}')

# Split interval equally using Python linspace
python3 -c "import numpy as np; [open(f'region_chunk_{i:02d}','w').write(f'$chrom\t{int(s)}\t{int(e)}\n') for i,(s,e) in enumerate(zip(np.linspace($first_pos,$last_pos,$CHUNKS+1)[:-1], np.linspace($first_pos,$last_pos,$CHUNKS+1)[1:]))]"

echo "Processing chunks in parallel..."
ls region_chunk_* | parallel -j $CHUNKS 'bcftools view '"$input_file"' -R {} -Ou | bcftools view -e '"'"'FILTER~"NO_HQ_GENOTYPES"'"'"' -Ou | bcftools annotate --set-id +'"'"'%CHROM\_%POS\_%REF\_%ALT'"'"' -Oz -o chunk{}.vcf.gz && echo {}'

echo "Concatenating chunks..."
bcftools concat -n -f <(ls chunk*.vcf.gz) -Oz -o tmp.vcf.gz && rm chunk*.vcf.gz region_chunk_*
tabix tmp.vcf.gz

echo "Done! Output file: tmp.vcf.gz"
echo "Original #variants"
file="$input_file"; bcftools index -s "$file" | awk '{sum+=$3} END{print sum}'
echo "Filtered #variants (after removing NO_HQ_GENOTYPES)"
file=tmp.vcf.gz; bcftools index -s "$file" | awk '{sum+=$3} END{print sum}'



