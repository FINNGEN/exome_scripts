#!/usr/bin/env python3
"""localize/bge.py — build the `gcloud storage cp` commands that copy a
completed bge.wdl (bge_qc) run's outputs into release/ (per
release/finngen_R14_exome_readme's "File structure" section):

  merged_vcf / merged_vcf_tbi -> data/qc_vcf_full/
  report                      -> documentation/

Per-chromosome intermediates (original_stats, filtered_stats,
validation_reports) aren't part of the release structure and are skipped —
only the final merged VCF+report matter.

Usage:
  python3 scripts/localize/bge.py <metadata.json> <gs://bucket/.../release/>
"""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import common

# (output_suffix, dest_subpath, fixed_filename_or_None, nested)
MAPPING = [
    ("merged_vcf",     "data/qc_vcf_full", None, False),
    ("merged_vcf_tbi", "data/qc_vcf_full", None, False),
    ("report",         "documentation",     None, False),
]

if __name__ == "__main__":
    common.main(MAPPING, __doc__)
