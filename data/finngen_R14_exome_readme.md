# FINNGEN EXOME DATA

## finngen_R14 exome data

Exome sequencing data processed from 45,399 samples across four sequencing batches, of which 43,302 were successfully mapped to existing FinnGen IDs. These data are personal data and must be treated according to the Finnish Personal Data Act 523/1999, EU Data Protection Directive 95/46/EC and EU General Data Protection Regulation (GDPR).

> For detailed pipeline documentation and source code, see the [GitHub repository](https://github.com/piotor87/exome_scripts).

The pipeline processes exome datasets from multiple sequencing batches through three main steps: quality control filtering, genetic verification of individual identity and ID mapping to existing FinnGen IDs, and LD computation against the FinnGen imputed genotype array.

| Dataset | Samples |
|---|---:|
| finngen_wes_gnomad_v4 | 25,201 |
| fimm-daly_finnish_gvs_bge_callset_1_padded_split_FINBBonly | 12,405 |
| THLBB2023_14_WES_Botnia | 7,164 |
| likely_pathogenic_annot_annotated_full_header_fix_resampled | 629 |
| **Total** | **45,399** |

---

### Processing and minimal QC

Each input exome dataset was independently processed through a per-chromosome QC pipeline run in parallel. The following operations were applied to each chromosome:

- **Chromosome name normalisation**: non-`chr`-prefixed contig names are renamed to the standard `chr` prefix.
- **FASTA normalisation**: variants are normalised against the GRCh38 reference using `bcftools norm`. Multi-allelic sites are split into biallelic records, indels are left-aligned, and REF mismatches are flagged and excluded.
- **Genotype masking**: genotypes with DP < 10 or GQ < 20 are set to missing.
- **Tag recalculation**: AC, AN and AF tags are recomputed after masking to reflect the updated genotype counts.
- **Variant filtering**: variants with AC = 0 after masking or spanning deletions (ALT="*") are removed.
- **Variant ID standardisation**: variant IDs are set to `CHROM_POS_REF_ALT` format.

Chromosomes are concatenated in their original order to produce a single QC-annotated VCF per dataset. A per-dataset report summarises the number of variants before and after filtering, by chromosome.

### Sample Renaming and ID Mapping

Exome sample IDs were mapped to FinnGen IDs through a genotype-based identity matching procedure. HM3 variants present in both the exome datasets and the FinnGen imputed array were used to compute pairwise KING kinship scores between all exome and FinnGen array samples. Each exome sample was assigned a list of candidate FinnGen IDs based on kinship, and resolution was applied as follows:

- **`ID_CONFIRMED`**: a single candidate was found and sample IDs match — unambiguous mapping.
- **`UNIQUE`**: a single candidate was found by genetics only, with no ID match.
- **`RESOLVED_BY_ID`**: multiple candidates found (e.g. twins in the reference); the query ID matched exactly one candidate.
- **`RESOLVED_BY_ALIAS`**: multiple candidates found, but all belong to the same alias group. Aliases are known FinnGen ID equivalences (e.g. individuals enrolled under different IDs across biobanks or sequencing batches) provided as an external file. The ambiguity is resolved by treating the group as a single identity.
- **`AMBIGUOUS_UNRESOLVED`**: multiple candidates with no alias resolution — sample excluded.
- **`MISSING`**: no genetic match found — sample excluded.
- **`HET_EXCLUDED`**: sample removed before KING due to high heterozygosity (F > 0.3); genotype data deemed unreliable for identity matching — sample excluded.
- **`CONFLICT_KEPT/DROPPED`**: when multiple exome samples from different datasets resolve to the same FinnGen ID, a priority tiebreak is applied (`ID_CONFIRMED` > `RESOLVED_BY_ALIAS` > `RESOLVED_BY_ID` > `UNIQUE`). The highest-priority match is kept; the rest are dropped.

The resulting mapping is stored in `finngen_R14_exome_id_mapping.tsv` (see File structure below). Each QC-annotated VCF is then subset to confirmed samples and reheadered with FinnGen IDs, producing a per-chromosome VCF per dataset. A flowchart visualising the resolution process across all datasets is available in the Documentation/Figures section.

#### ID mapping summary

| Group | Status | N | Description |
|---|---|---:|---|
| MATCHED | `ID_CONFIRMED` | 40,189 | Single candidate; confirmed by matching IDs |
| MATCHED | `RESOLVED_BY_ALIAS` | 1,524 | Multiple candidates resolved via known alias group |
| MATCHED | `CONFLICT_KEPT` (was `ID_CONFIRMED`) | 1,347 | Contested FinnGen ID; won priority tiebreak |
| MATCHED | `RESOLVED_BY_ID` | 160 | Twins in ref; query ID matched one candidate |
| MATCHED | `UNIQUE` | 52 | Single candidate; matched by genetics only |
| MATCHED | `CONFLICT_KEPT` (was `RESOLVED_BY_ALIAS`) | 30 | Contested FinnGen ID; won priority tiebreak |
| **Total matched** | | **43,302** | |
| DROPPED | `CONFLICT_DROPPED` (was `ID_CONFIRMED`) | 1,146 | Contested FinnGen ID; lost tiebreak |
| DROPPED | `CONFLICT_DROPPED` (was `RESOLVED_BY_ALIAS`) | 286 | Contested FinnGen ID; lost tiebreak |
| DROPPED | `CONFLICT_DROPPED` (was `UNIQUE`) | 18 | Contested FinnGen ID; lost tiebreak |
| DROPPED | `CONFLICT_DROPPED` (was `RESOLVED_BY_ID`) | 6 | Contested FinnGen ID; lost tiebreak |
| **Total dropped** | | **1,456** | |
| NO MATCH | `MISSING` | 636 | No genetic match found |
| EXCLUDED | `HET_EXCLUDED` | 5 | Excluded by heterozygosity filter (F > 0.3) prior to KING — WES dataset only |
| **Total** | | **45,399** | |

### LD

LD (r²) between FinnGen imputed array SNPs and exome variants is computed genome-wide. FinnGen and exome data are first converted to plink2 format and merged per chromosome. plink2 `--r2-unphased` is then used to compute LD between each FinnGen SNP and all exome variants within a fixed genomic window. Each SNP pair is annotated with coding/non-coding status for both variants. Per-chromosome results are concatenated into a single genome-wide file. A filtered version retaining only pairs with r² ≥ [min_r2] is produced alongside per-chromosome summary statistics and figures.

---

## File structure

| File | Description |
|---|---|
| `data/finngen_R14_exome_id_mapping.tsv` | QRY→REF sample ID mapping with resolution status for all datasets |
| `data/qc_vcf_full/[EXOME_DATASET].QC_ANNOTATED.vcf.gz` | QC-filtered and normalised VCF for each input dataset |
| `data/qc_vcf_full/[EXOME_DATASET].QC_ANNOTATED.vcf.gz.tbi` | Index for QC-filtered VCF |
| `data/renamed_vcf_chr/[EXOME_DATASET].QC_ANNOTATED_fg_ids_chr[N].vcf.gz` | Per-chromosome VCF with FinnGen sample IDs |
| `data/renamed_vcf_chr/[EXOME_DATASET].QC_ANNOTATED_fg_ids_chr[N].vcf.gz.tbi` | Index for renamed VCF |
| `data/finngen_R14_exome.ld.tsv.gz` | All genome-wide FinnGen–exome LD pairs with coding annotation |
| `data/finngen_R14_exome_r[min_r2]_ld.tsv.gz` | LD pairs filtered to r² ≥ [min_r2] |

### Logs

| File | Description |
|---|---|
| `documentation/[EXOME_DATASET].QC_ANNOTATED.report.txt` | Per-dataset QC filtering statistics by chromosome |
| `documentation/finngen_R14_exome_r[min_r2]_ld_stats.tsv` | Per-chromosome LD summary statistics |

### Figures

| File | Description |
|---|---|
| `documentation/finngen_R14_exome_id_mapping_flowchart.png` | Flowchart of sample ID resolution across datasets |
| `documentation/finngen_R14_exome_r[min_r2]_fig1_variants.png` | Unique variants per chromosome, stacked coding/non-coding |
| `documentation/finngen_R14_exome_r[min_r2]_fig2_pairs.png` | LD pair breakdown per chromosome |
| `documentation/finngen_R14_exome_r[min_r2]_fig3_r2_dist.png` | r² distribution by coding category |
| `documentation/finngen_R14_exome_r[min_r2]_fig4_coding_frac.png` | Coding fraction per chromosome |
