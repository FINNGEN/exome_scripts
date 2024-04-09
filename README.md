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

First one needs to run `exome_convert.wdl` to get a decent subet after filtering for DP/GQ and annotating to CHROM_POS_REF_ALT

Then one needs to convert to plink to shared FG variants and run king to tget a `kin` and `kin0` file. One needs to reminds also to convert the plink files to `--maj-ref` with plink1.


The one can run `assign_ids.wdl` to assign the mapping from/to exome/fg IDs.

Then one can run again the exome_convert.wdl with different json to make sure to subset and reheader the samples