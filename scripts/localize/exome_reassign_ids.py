#!/usr/bin/env python3
"""localize/exome_reassign_ids.py — build the `gcloud storage cp` commands
that copy a completed exome_reassign_ids.wdl run's outputs into release/
(per release/finngen_R14_exome_readme's "File structure" section):

  id_mapping                             -> data/finngen_R14_exome_id_mapping.tsv (fixed name, top of data/)
  chrom_vcfs[] / chrom_tbis[]            -> data/renamed_vcf_chr/
  hwe_beds[] / hwe_bims[] / hwe_fams[]   -> data/renamed_plink_chr/
  dataset_hwe_summaries[]                -> documentation/hwe/

Rename-stage and HWE-stage outputs are only present when run_rename = true
and are skipped (with a note on stderr) otherwise.

combined_summary / combined_plot / id_mapping_stats / id_mapping_md /
id_mapping_flowchart are NOT copied here — those go to this repo's own
data/ folder for documentation, not to the release bucket.

Usage:
  python3 scripts/localize/exome_reassign_ids.py <metadata.json> <gs://bucket/.../release/>
"""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import common

# (output_suffix, dest_subpath, fixed_filename_or_None, nested)
MAPPING = [
    ("id_mapping",            "data", "finngen_R14_exome_id_mapping.tsv", False),
    ("chrom_vcfs",             "data/renamed_vcf_chr",   None, False),
    ("chrom_tbis",             "data/renamed_vcf_chr",   None, False),
    ("hwe_beds",               "data/renamed_plink_chr", None, False),
    ("hwe_bims",               "data/renamed_plink_chr", None, False),
    ("hwe_fams",               "data/renamed_plink_chr", None, False),
    ("dataset_hwe_summaries",  "documentation/hwe",      None, False),
]

if __name__ == "__main__":
    common.main(MAPPING, __doc__)
