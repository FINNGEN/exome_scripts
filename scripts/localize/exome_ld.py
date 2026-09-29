#!/usr/bin/env python3
"""localize/exome_ld.py — build the `gcloud storage cp` commands that copy a
completed exome_ld.wdl run's outputs into release/ (per
release/finngen_R14_exome_readme's "File structure" section):

  merged_plink (bed/bim/fam) + afreq -> data/plink_fg_merged_chr/
  annotated_ld -> data/finngen_R14_exome.ld_annotated.tsv.gz  (fixed name, top of data/)
  raw_ld       -> data/finngen_R14_exome.ld.tsv.gz             (fixed name, top of data/)

merged_plink is Array[Array[File]] (one [bed,bim,fam] triple per chrom) —
flattened to a single file list before copying. Per-chrom intermediates
(ld_results, merged_vcf) and leads_tsv_out aren't part of the release
structure and are skipped.

Usage:
  python3 scripts/localize/exome_ld.py <metadata.json> <gs://bucket/.../release/>
"""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import common

# (output_suffix, dest_subpath, fixed_filename_or_None, nested)
MAPPING = [
    ("merged_plink", "data/plink_fg_merged_chr", None, True),   # Array[Array[File]] -> flatten
    ("afreq",        "data/plink_fg_merged_chr", None, False),
    ("annotated_ld", "data", "finngen_R14_exome.ld_annotated.tsv.gz", False),
    ("raw_ld",       "data", "finngen_R14_exome.ld.tsv.gz",           False),
]

if __name__ == "__main__":
    common.main(MAPPING, __doc__)
