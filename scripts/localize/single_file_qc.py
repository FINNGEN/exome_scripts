#!/usr/bin/env python3
"""localize/single_file_qc.py — build the `gcloud storage cp` commands that
copy a completed single_file_qc.wdl run's outputs into release/ (per
release/finngen_R14_exome_readme's "File structure" section):

  filtered_vcfs[] / filtered_vcf_tbis[] -> data/qc_vcf_full/
  validation_reports[]                  -> documentation/

Unlike gnomad_wes_finns.wdl/bge.wdl, this workflow has no merge step — its
outputs stay per-chromosome, so these are arrays copied as-is (basenames
preserved).

Usage:
  python3 scripts/localize/single_file_qc.py <metadata.json> <gs://bucket/.../release/>
"""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import common

# (output_suffix, dest_subpath, fixed_filename_or_None, nested)
MAPPING = [
    ("filtered_vcfs",      "data/qc_vcf_full", None, False),
    ("filtered_vcf_tbis",  "data/qc_vcf_full", None, False),
    ("validation_reports", "documentation",     None, False),
]

if __name__ == "__main__":
    common.main(MAPPING, __doc__)
