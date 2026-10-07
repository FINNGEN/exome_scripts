#!/bin/bash
set -o xtrace


DATA_DIR='/home/pete/fg-3/exome/renamed_final/vcf/'
OUT_DIR='/mnt/disks/data/exome/qc/'

echo $DATA_DIR $OUT_DIR
bcftools query -v <(ls $DATA_DIR/*gz |sort -V)   -f   '[%CHROM\t%POS\t%SAMPLE\tAD=%AD,DP=%DP,GT=%GT\n]' -i '(GT="het" & DP<20 & AD[:1]/DP>.8) | (GT="het" & DP<20 & AD[:1]/DP<.2)  ' | bgzip -c > $OUT_DIR"AB.gz" &
bcftools query -v <(ls $DATA_DIR/*gz |sort -V)   -f   '[%CHROM\t%POS\t%SAMPLE\tAD=%AD,DP=%DP,GT=%GT\n]' -i '(GT="het" & DP<20)' | bgzip -c > $OUT_DIR"DP.gz"


	      
    
    
	 
