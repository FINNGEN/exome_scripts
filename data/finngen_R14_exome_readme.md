# FINNGEN EXOME DATA

## finngen_R14 exome data

Exome sequencing data processed from 45,399 samples across four sequencing batches, of which 43,302 were successfully mapped to existing FinnGen IDs. These data are personal data and must be treated according to the Finnish Personal Data Act 523/1999, EU Data Protection Directive 95/46/EC and EU General Data Protection Regulation (GDPR).

> For detailed pipeline documentation and source code, see the [GitHub repository](https://github.com/FINNGEN/exome_scripts)

The pipeline processes exome datasets from multiple sequencing batches through three main steps: quality control filtering, genetic verification of individual identity and ID mapping to existing FinnGen IDs, and merging with the FinnGen imputed genotype array into a joint plink dataset.

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

### Annotated LD

For each FinnGen credible set lead variant, exome variants in LD (r² ≥ 0.05, within a 1 Mb window) are reported alongside phenotype and credible set annotation. The file has one row per (FG lead variant, exome variant, phenotype) triple — the same LD pair is duplicated across all phenotypes for which that FG variant is a credible set lead. Columns:

| Column | Description |
|---|---|
| `FG_SNP` | FinnGen array variant (credible set lead) |
| `EXOME_SNP` | Exome variant in LD with the lead |
| `R2` | Unphased r² between the pair |
| `exome_consequence` | Most severe VEP consequence for the exome variant (`NA` if non-coding) |
| `EXOME_AF` | Allele frequency of the exome variant in the merged dataset |
| `PHENO` | FinnGen phenotype abbreviation |
| `lead_mlogp` | –log₁₀(p) of the lead variant for this phenotype |
| `lead_beta` | Effect size (beta) of the lead variant |
| `lead_af_alt` | Allele frequency of the lead variant in FinnGen |
| `good_cs` | Whether the credible set passes quality criteria |
| `cs_type` | Credible set classification: `coding` (lead has a coding variant), `functional_relaxed` (relaxed functional variant but no coding), or `NA` (no functional explanation) |
| `functional_var` | Top functional variant in the credible set (first entry of `functional_variants_relaxed`); `NA` if none |
| `functional_var_r2` | r² between `functional_var` and the lead; `NA` if none |

### Merged plink data

FinnGen imputed array variants and exome variants are jointly converted to plink1 binary format and merged per chromosome. For each chromosome, the FG plink fileset (subsetted to the 43,302 matched exome samples) and all four exome plink filesets are merged via `plink --merge-list` into a single combined BED/BIM/FAM. Exome-private variants are included alongside FinnGen array variants; variants already present in the FG BIM are excluded from the exome filesets to avoid duplication. The resulting filesets span chromosomes 1–23 and contain ~31.1 million variants across 43,302 samples.

---

## File structure

### Data

| File | Description |
|---|---|
| `finngen_R14_exome_id_mapping.tsv` | QRY→REF sample ID mapping with resolution status for all datasets |
| `qc_vcf_full/[EXOME_DATASET].QC_ANNOTATED.vcf.gz` | QC-filtered and normalised VCF for each input dataset |
| `qc_vcf_full/[EXOME_DATASET].QC_ANNOTATED.vcf.gz.tbi` | Index for QC-filtered VCF |
| `renamed_vcf_chr/[EXOME_DATASET].QC_ANNOTATED_fg_ids_chr[N].vcf.gz` | Per-chromosome VCF with FinnGen sample IDs |
| `renamed_vcf_chr/[EXOME_DATASET].QC_ANNOTATED_fg_ids_chr[N].vcf.gz.tbi` | Index for renamed VCF |
| `finngen_R14_exome.ld_annotated.tsv.gz` | Exome variants in LD with FinnGen credible set leads, annotated with phenotype and credible set metadata (~900k rows) |
| `plink_fg_merged_chr/finngen_R14_exome_chr[N].bed` | Per-chromosome merged plink BED (FG array + exome variants, 43,302 samples, chr1–23) |
| `plink_fg_merged_chr/finngen_R14_exome_chr[N].bim` | BIM file; variant IDs in `CHROM_POS_REF_ALT` format; ~31.1 M variants genome-wide |
| `plink_fg_merged_chr/finngen_R14_exome_chr[N].fam` | FAM file with FinnGen sample IDs |

### Documentation

#### Logs

| File | Description |
|---|---|
| `[EXOME_DATASET].QC_ANNOTATED.report.txt` | Per-dataset QC filtering statistics by chromosome |

#### Figures

| File | Description |
|---|---|
| `FG_EXOME_resolved_flowchart.png` | Flowchart of sample ID resolution across datasets |
