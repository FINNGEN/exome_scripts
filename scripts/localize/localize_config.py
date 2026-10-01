"""localize_config.py — per-workflow output -> release/ mapping tables.

Each WDL's outputs map onto the release/ file structure documented in
release/finngen_R14_exome_readme ("## File structure" section is the ground
truth for every path used here). Keyed by workflowName — the `workflow
<name> { ... }` declared inside each WDL, NOT the filename (e.g.
wdl/gnomad_wes_finns.wdl declares `workflow gnomad_wes_finns_chrom`) — so
localize.py can auto-detect which mapping to use from a run's metadata JSON.

Each mapping is a list of (output_suffix, dest_subpath, fixed_filename_or_None, nested):
  output_suffix:  the WDL output name, matched by suffix against the
                  workflow-name-prefixed keys in metadata's "outputs"
                  (e.g. "chrom_vcfs" matches "exome_reassign_ids.chrom_vcfs")
  dest_subpath:   path under release/ to copy into
  fixed_filename_or_None: for a scalar output renamed to a fixed name on
                  copy (e.g. id_mapping -> finngen_R14_exome_id_mapping.tsv);
                  None means an array output where each file's own basename
                  is kept
  nested:         True if the output is Array[Array[File]] and needs
                  flattening first (only exome_ld's merged_plink today)
"""

WORKFLOWS = {

    # wdl/gnomad_wes_finns.wdl
    "gnomad_wes_finns_chrom": [
        ("concatenated_vcf",     "data/qc_vcf_full", None, False),
        ("concatenated_vcf_tbi", "data/qc_vcf_full", None, False),
        ("report",               "documentation",     None, False),
        # Per-chromosome intermediates (filtered_vcfs, filtered_vcf_tbis,
        # original_stats, filtered_stats, validation_reports) aren't part of
        # the release structure and are intentionally omitted.
    ],

    # wdl/bge.wdl (workflow bge_qc)
    "bge_qc": [
        ("merged_vcf",     "data/qc_vcf_full", None, False),
        ("merged_vcf_tbi", "data/qc_vcf_full", None, False),
        ("report",         "documentation",     None, False),
        # original_stats/filtered_stats/validation_reports: intermediates, omitted.
    ],

    # wdl/single_file_qc.wdl — no merge step, outputs stay per-chromosome
    "single_file_qc": [
        ("filtered_vcfs",      "data/qc_vcf_full", None, False),
        ("filtered_vcf_tbis",  "data/qc_vcf_full", None, False),
        ("validation_reports", "documentation",     None, False),
    ],

    # wdl/exome_ld.wdl
    "exome_ld": [
        ("merged_plink", "data/plink_fg_merged_chr", None, True),   # Array[Array[File]] -> flatten
        ("afreq",        "data/plink_fg_merged_chr", None, False),
        ("annotated_ld", "data", "finngen_R14_exome.ld_annotated.tsv.gz", False),
        ("raw_ld",       "data", "finngen_R14_exome.ld.tsv.gz",           False),
        # ld_results (per-chrom), merged_vcf, leads_tsv_out: intermediates, omitted.
    ],

    # wdl/exome_reassign_ids.wdl
    "exome_reassign_ids": [
        ("id_mapping",            "data",                   "finngen_R14_exome_id_mapping.tsv", False),
        ("chrom_vcfs",             "data/renamed_vcf_chr",   None, False),
        ("chrom_tbis",             "data/renamed_vcf_chr",   None, False),
        ("hwe_beds",               "data/renamed_plink_chr", None, False),
        ("hwe_bims",               "data/renamed_plink_chr", None, False),
        ("hwe_fams",               "data/renamed_plink_chr", None, False),
        ("dataset_hwe_summaries",  "documentation",          None, False),
        # combined_summary/combined_plot/id_mapping_stats/id_mapping_md/
        # id_mapping_flowchart go to this repo's own data/ folder for
        # documentation, not to the release bucket — intentionally omitted.
        # Only present when run_rename = true; missing otherwise (skipped
        # with a note on stderr, not an error).
    ],

}
