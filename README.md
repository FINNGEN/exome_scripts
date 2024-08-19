# EXOME data processing

## kinship.sh

Runs through output of FG vs EXOME data and it procudes summaries of correct/incorrect ID mappings

Usage:
```
bash kinship.sh ./mapping_test/ mapping.kin
```

Check for paths to kinship data in the script. 

## assign_id.py
Produces a file that can be passed to plink `--update-ids` in order to update the exome sample ids

Usage:

```
python3 ~/Dropbox/Projects/misc/exome/assign_id.py /mnt/disks/data/exome/kinship/full_exome/ /mnt/disks/data/exome/plink/exome_merged.fam
```


## WDL

First one needs to run `exome_annotate.wdl` to get a decent subet after filtering for DP/GQ and annotating to CHROM_POS_REF_ALT

Then one needs to convert to plink to shared FG variants and run king to get a `kin` and `kin0` file. One needs to reminds also to convert the plink files to `--maj-ref` with plink1.

The one can run `assign_ids.wdl` to assign the mapping from/to exome/fg IDs.

After that one can run 'exome_convert.wdl` to build the vcf/plink/bgen files with the new sample ids

`exome_ld.wdl` is the wdl where plink FG data is subset to exome samples and exome data is subset to only non FG snps, so that they can be merged and ld run.

