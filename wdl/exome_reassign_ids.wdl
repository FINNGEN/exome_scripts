version 1.0

# Chains exome_duplicates.wdl (kinship-based ID resolution) and
# exome_rename.wdl (subset + rename to FinnGen IDs) into one submission.
#
# Run with run_rename = false first to compute and inspect the ID mapping
# (id_mapping / id_mapping_stats / id_mapping_flowchart) without paying the
# cost of the rename step. Once the mapping looks right, resubmit the same
# workflow with only run_rename = true changed — Cromwell call-caches every
# upstream call (none of their inputs changed), so only the newly enabled
# exome_rename sub-workflow actually runs.

import "exome_duplicates.wdl" as dup_wf
import "exome_rename.wdl" as rename_wf

workflow exome_reassign_ids {
  input {
    # ---- shared: same [DATASET, QC_ANNOTATED.vcf.gz] pairs feed both stages ----
    Array[Array[String]] vcf_pairs

    # ---- exome_duplicates (kinship / ID resolution) ----
    File   plink_bed
    String plink_prefix
    Int    n_regions     = 100
    File?  berisa_blocks
    File?  aliases
    Int    target_snps   = 10000
    Float  max_het_F     = 0.3
    Int    chunk_size    = 10000
    String out_prefix    = "finngen_R14_exome"

    # ---- exome_rename (subset + rename to FinnGen IDs) ----
    Boolean       run_rename = false   # flip to true and resubmit once id_mapping looks right
    Array[String] chroms = [
      "chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
      "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
      "chr20","chr21","chr22","chrX","chrY"
    ]
    Int    chunk_mb = 500
    String suffix   = "fg_ids"
  }

  call dup_wf.exome_duplicates {
    input:
      vcf_pairs     = vcf_pairs,
      plink_bed     = plink_bed,
      plink_prefix  = plink_prefix,
      n_regions     = n_regions,
      berisa_blocks = berisa_blocks,
      aliases       = aliases,
      target_snps   = target_snps,
      max_het_F     = max_het_F,
      chunk_size    = chunk_size,
      out_prefix    = out_prefix
  }

  if (run_rename) {
    call rename_wf.exome_rename {
      input:
        resolved_mapping = exome_duplicates.id_mapping,
        vcf_pairs        = vcf_pairs,
        chroms           = chroms,
        chunk_mb         = chunk_mb,
        suffix           = suffix
    }
  }

  output {
    # ---- from exome_duplicates — always produced ----
    File combined_summary     = exome_duplicates.combined_summary
    File combined_plot        = exome_duplicates.combined_plot
    File id_mapping           = exome_duplicates.id_mapping
    File id_mapping_stats     = exome_duplicates.id_mapping_stats
    File id_mapping_md        = exome_duplicates.id_mapping_md
    File id_mapping_flowchart = exome_duplicates.id_mapping_flowchart

    # ---- from exome_rename — only present when run_rename = true ----
    Array[File]? chrom_vcfs = exome_rename.chrom_vcfs
    Array[File]? chrom_tbis = exome_rename.chrom_tbis
  }
}
