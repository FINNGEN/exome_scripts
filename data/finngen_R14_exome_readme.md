# FINNGEN EXOME DATA

## finngen_R14 exome data

Exome sequencing data for the Analysis Team and collaborating partners. These data are personal data and must be treated according to the Finnish Personal Data Act 523/1999, EU Data Protection Directive 95/46/EC and EU General Data Protection Regulation (GDPR).

> For detailed pipeline documentation and source code, see the [GitHub repository](https://github.com/piotor87/exome_scripts).

The pipeline processes exome datasets from multiple sequencing batches through three main steps: quality control filtering, sample renaming to FinnGen IDs, and LD computation against the FinnGen imputed genotype array.

---

### Quality Control

Each input exome dataset was independently processed through a per-chromosome QC pipeline run in parallel. The following operations were applied to each chromosome:

- **Chromosome name normalisation**: non-`chr`-prefixed contig names are renamed to the standard `chr` prefix.
- **FASTA normalisation**: variants are normalised against the GRCh38 reference using `bcftools norm`. Multi-allelic sites are split into biallelic records, indels are left-aligned, and REF mismatches are flagged and excluded.
- **Genotype masking**: genotypes failing per-call quality criteria are set to missing.
- **Tag recalculation**: AC, AN and AF tags are recomputed after masking to reflect the updated genotype counts.
- **Variant filtering**: variants failing quality thresholds (e.g. AC=0 after masking) are removed.
- **Variant ID standardisation**: variant IDs are set to `CHROM_POS_REF_ALT` format.

Chromosomes are concatenated in their original order to produce a single QC-annotated VCF per dataset. A per-dataset report summarises the number of variants before and after filtering, by chromosome.

### Sample Renaming and ID Mapping

Exome sample IDs were mapped to FinnGen IDs through a genotype-based identity matching procedure. HM3 variants present in both the exome datasets and the FinnGen imputed array were used to compute pairwise KING kinship scores between all exome and FinnGen array samples. Each exome sample was assigned a list of candidate FinnGen IDs based on kinship, and resolution was applied as follows:

- **`ID_CONFIRMED`**: a single candidate was found and sample IDs match — unambiguous mapping.
- **`UNIQUE`**: a single candidate was found by genetics only, with no ID match.
- **`RESOLVED_BY_ID`**: multiple candidates found (e.g. twins in the reference); the query ID matched exactly one candidate.
- **`RESOLVED_BY_ALIAS`**: multiple candidates found, but all belong to the same alias group. Aliases are known FinnGen ID equivalences (e.g. individuals enrolled under different IDs across biobanks or sequencing batches) provided as an external file. The ambiguity is resolved by treating the group as a single identity.
- **`AMBIGUOUS_UNRESOLVED`**: multiple candidates with no alias resolution — sample excluded.
- **`MISSING`**: no kinship match found — sample excluded.
- **`CONFLICT_KEPT/DROPPED`**: when multiple exome samples from different datasets resolve to the same FinnGen ID, a priority tiebreak is applied (`ID_CONFIRMED` > `RESOLVED_BY_ALIAS` > `RESOLVED_BY_ID` > `UNIQUE`). The highest-priority match is kept; the rest are dropped.

The resulting mapping is stored in `/home/pete/fg-3/exome_v2/release/data/finngen_R14_exome_id_mapping.tsv`. Each QC-annotated VCF is then subset to confirmed samples and reheadered with FinnGen IDs, producing a per-chromosome VCF per dataset. A flowchart visualising the resolution process across all datasets is available in the Documentation/Figures section.

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
| NO MATCH | `MISSING` | 636 | No KING match found |
| **Total** | | **45,394** | |

### LD

LD (r²) between FinnGen imputed array SNPs and exome variants is computed genome-wide. FinnGen and exome data are first converted to plink2 format and merged per chromosome. plink2 `--r2-unphased` is then used to compute LD between each FinnGen SNP and all exome variants within a fixed genomic window. Each SNP pair is annotated with coding/non-coding status for both variants. Per-chromosome results are concatenated into a single genome-wide file. A filtered version retaining only pairs with r² ≥ [min_r2] is produced alongside per-chromosome summary statistics and figures.

---

## File structure

### Data

| File | Description |
|---|---|
| `/home/pete/fg-3/exome_v2/release/data/finngen_R14_exome_id_mapping.tsv` | QRY→REF sample ID mapping with resolution status for all datasets |
| `qc_vcf_full/[EXOME_DATASET].QC_ANNOTATED.vcf.gz` | QC-filtered and normalised VCF for each input dataset |
| `qc_vcf_full/[EXOME_DATASET].QC_ANNOTATED.vcf.gz.tbi` | Index for QC-filtered VCF |
| `renamed_vcf_chr/[EXOME_DATASET].QC_ANNOTATED_fg_ids_chr[N].vcf.gz` | Per-chromosome VCF with FinnGen sample IDs |
| `renamed_vcf_chr/[EXOME_DATASET].QC_ANNOTATED_fg_ids_chr[N].vcf.gz.tbi` | Index for renamed VCF |
| `finngen_R14_exome.ld.tsv.gz` | All genome-wide FinnGen–exome LD pairs with coding annotation |
| `finngen_R14_exome_r[min_r2]_ld.tsv.gz` | LD pairs filtered to r² ≥ [min_r2] |
| `finngen_R14_exome_r[min_r2]_ld_stats.tsv` | Per-chromosome LD summary statistics |

### Documentation

#### Logs

| File | Description |
|---|---|
| `[EXOME_DATASET].QC_ANNOTATED.report.txt` | Per-dataset QC filtering statistics by chromosome |

#### Figures

| File | Description |
|---|---|
| `FG_EXOME_resolved_flowchart.png` | Flowchart of sample ID resolution across datasets |
| `finngen_R14_exome_r[min_r2]_fig1_variants.png` | Unique variants per chromosome, stacked coding/non-coding |
| `finngen_R14_exome_r[min_r2]_fig2_pairs.png` | LD pair breakdown per chromosome |
| `finngen_R14_exome_r[min_r2]_fig3_r2_dist.png` | r² distribution by coding category |
| `finngen_R14_exome_r[min_r2]_fig4_coding_frac.png` | Coding fraction per chromosome |
