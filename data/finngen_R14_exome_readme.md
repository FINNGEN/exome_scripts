# FINNGEN EXOME DATA

## finngen_R14 exome data

Exome sequencing data processed from 45,375 samples across four sequencing batches, of which 43,289 were successfully mapped to existing FinnGen IDs. These data are personal data and must be treated according to the Finnish Personal Data Act 523/1999, EU Data Protection Directive 95/46/EC and EU General Data Protection Regulation (GDPR).

> For detailed pipeline documentation and source code, see the [GitHub repository](https://github.com/FINNGEN/exome_scripts)

The pipeline processes exome datasets from multiple sequencing batches through three main steps: quality control filtering, genetic verification of individual identity and ID mapping to existing FinnGen IDs, and merging with the FinnGen imputed genotype array into a joint plink dataset.
The exomes do have 1-on-1 mapping with imputed FinnGen participants with the same FinnGen ID

| Long name | Short name | Samples | Description |
|---|---|---:|---|
| gnomAD v4 Finns subset | `gnomad_wes_finns` | 25,197 | GnomAD Finns subset of individuals  already in FinnGen |
| Blended Genome Exome scizophrenia, bipolar, controls | `BGE_scz_bp_ctrl` | 12,389 | Broad Institute sequenced Bipolar & schizophrenia and shared control cohort |
| Botnia THL diabetes study | `botnia` | 7,161 | Botnia Diabetes cohort |
| Autosomal dominant polycystic kidney disease WES (ADPKD) | `ADPKD` | 628 | Autosomal dominant polycystic kidney disease patients |
| **Total** | | **45,375** | |

---

### Processing and minimal QC

Each input exome dataset was independently processed through a per-chromosome QC pipeline run in parallel. The following operations were applied to each chromosome:

- **Sample exclusion**: samples on a registry-mandated denial list — expanded to include every known alias of each denied ID — are removed before any other QC step. 24 samples were excluded this way across the four datasets (BGE_scz_bp_ctrl 16, gnomad_wes_finns 4, botnia 3, ADPKD 1); the sample counts in the table above are post-exclusion.
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
| MATCHED | `ID_CONFIRMED` | 40,177 | Single candidate; confirmed by matching IDs |
| MATCHED | `RESOLVED_BY_ALIAS` | 1,523 | Multiple candidates resolved via known alias group |
| MATCHED | `CONFLICT_KEPT` (was `ID_CONFIRMED`) | 1,347 | Contested FinnGen ID; won priority tiebreak |
| MATCHED | `RESOLVED_BY_ID` | 160 | Twins in ref; query ID matched one candidate |
| MATCHED | `UNIQUE` | 52 | Single candidate; matched by genetics only |
| MATCHED | `CONFLICT_KEPT` (was `RESOLVED_BY_ALIAS`) | 30 | Contested FinnGen ID; won priority tiebreak |
| **Total matched** | | **43,289** | |
| DROPPED | `CONFLICT_DROPPED` (was `ID_CONFIRMED`) | 1,146 | Contested FinnGen ID; lost tiebreak |
| DROPPED | `CONFLICT_DROPPED` (was `RESOLVED_BY_ALIAS`) | 286 | Contested FinnGen ID; lost tiebreak |
| DROPPED | `CONFLICT_DROPPED` (was `UNIQUE`) | 18 | Contested FinnGen ID; lost tiebreak |
| DROPPED | `CONFLICT_DROPPED` (was `RESOLVED_BY_ID`) | 6 | Contested FinnGen ID; lost tiebreak |
| **Total dropped** | | **1,456** | |
| NO MATCH | `MISSING` | 625 | No genetic match found |
| EXCLUDED | `HET_EXCLUDED` | 5 | Excluded by heterozygosity filter (F > 0.3) prior to KING — WES dataset only |
| **Total** | | **45,375** | |

### FG/Exome LD

LD between FinnGen array variants and exome-private variants is computed genome-wide (r² ≥ 0.05, within a 1 Mb window) using plink2 `--r2-unphased` on the merged FG + exome plink dataset. Two files are produced — a raw LD file and a version annotated with FinnGen credible set data — both covering all chromosomes 1–23 in a single file.

#### Raw LD (`finngen_R14_exome.ld.tsv.gz`)

All FG→exome LD pairs passing the r² threshold. Columns:

| Column | Description |
|---|---|
| `FG_SNP` | FinnGen array variant |
| `EXOME_SNP` | Exome variant in LD |
| `R2` | Unphased r² between the pair |
| `exome_consequence` | Most severe VEP consequence for the exome variant |
| `EXOME_AF` | Allele frequency of the exome variant in the merged dataset |
| `exome_nearest_gene` | Gene corresponding to the exome variant's most severe VEP consequence |

#### Annotated LD (`finngen_R14_exome.ld_annotated.tsv.gz`)

For each FinnGen credible set lead variant, exome variants in LD (r² ≥ 0.05, within a 1 Mb window) are reported alongside phenotype and credible set annotation. The file has one row per (FG lead variant, exome variant, phenotype) triple — the same LD pair is duplicated across all phenotypes for which that FG variant is a credible set lead. Columns:

| Column | Description |
|---|---|
| `FG_SNP` | FinnGen array variant (credible set lead) |
| `EXOME_SNP` | Exome variant in LD with the lead |
| `R2` | Unphased r² between the pair |
| `exome_consequence` | Most severe VEP consequence for the exome variant (`NA` if non-coding) |
| `EXOME_AF` | Allele frequency of the exome variant in the merged dataset |
| `exome_nearest_gene` | Gene corresponding to the exome variant's most severe VEP consequence |
| `PHENO` | FinnGen phenotype abbreviation |
| `lead_mlogp` | –log₁₀(p) of the lead variant for this phenotype |
| `lead_beta` | Effect size (beta) of the lead variant |
| `lead_af_alt` | Allele frequency of the lead variant in FinnGen |
| `good_cs` | Whether the credible set passes quality criteria |
| `cs_log_bayes_factor` | Log Bayes factor for including this credible set vs. excluding it; collapses (≪ 1) for non-independent CS artifacts |
| `cs_type` | Credible set classification: `coding` (lead has a coding variant), `functional_relaxed` (relaxed functional variant but no coding), or `NA` (no functional explanation) |
| `functional_var` | Top functional variant in the credible set (first entry of `functional_variants_relaxed`); `NA` if none |
| `functional_var_r2` | r² between `functional_var` and the lead; `NA` if none |

### Merged plink data

FinnGen imputed array variants and exome variants are jointly converted to plink1 binary format and merged per chromosome. For each chromosome, the FG plink fileset (subsetted to the 43,289 matched exome samples) and all four exome plink filesets are merged via `plink --merge-list` into a single combined BED/BIM/FAM. Exome-private variants are included alongside FinnGen array variants; variants already present in the FG BIM are excluded from the exome filesets to avoid duplication. The resulting filesets span chromosomes 1–23 and contain ~31.1 million variants across 43,289 samples.

### Hardy-Weinberg equilibrium

For each exome dataset, Hardy-Weinberg equilibrium is tested per autosomal chromosome (chr1–chr22) on that dataset's own per-chromosome plink fileset using plink2 `--hardy`. Allele frequency (`ALT_FREQS`/`OBS_CT`) is computed directly from `--hardy`'s own genotype counts rather than a separate `--freq` run, then concatenated across chromosomes in order into one gzipped per-dataset summary.

#### Per-dataset HWE summary (`[EXOME_DATASET].hwe_summary.tsv.gz`)

One row per variant, autosomes only. Columns:

| Column | Description |
|---|---|
| `#CHROM` | Chromosome |
| `ID` | Variant ID (`CHROM_POS_REF_ALT` format) |
| `A1` | Reference allele for the HWE test |
| `AX` | Non-A1 allele(s) |
| `HOM_A1_CT` | Count of samples homozygous for A1 |
| `HET_A1_CT` | Count of heterozygous samples |
| `TWO_AX_CT` | Count of samples homozygous for the non-A1 allele |
| `O(HET_A1)` | Observed heterozygous A1 frequency |
| `E(HET_A1)` | Expected heterozygous A1 frequency under HWE |
| `P` | Hardy-Weinberg exact test p-value |
| `ALT_FREQS` | Alternate allele frequency (derived from `--hardy`'s own genotype counts) |
| `OBS_CT` | Number of allele observations (derived from `--hardy`'s own genotype counts) |

---

## File structure

### Data

| File | Description |
|---|---|
| `finngen_R14_exome_id_mapping.tsv` | QRY→REF sample ID mapping with resolution status for all datasets |
| `qc_vcf_full/[EXOME_DATASET].QC_ANNOTATED.vcf.gz[.tbi]` | QC-filtered and normalised VCF for each input dataset |
| `renamed_vcf_chr/[EXOME_DATASET].QC_ANNOTATED_fg_ids_chr[N].vcf.gz[.tbi]` | Per-chromosome VCF with FinnGen sample IDs |
| `renamed_plink_chr/[EXOME_DATASET].QC_ANNOTATED_fg_ids_chr[N].{bed,bim,fam}` | Per-chromosome plink dataset, autosomes only (chr1–22); same naming as the matching renamed VCF |
| `finngen_R14_exome.ld.tsv.gz` | All FG→exome LD pairs (r² ≥ 0.05, 1 Mb window), genome-wide; annotated with VEP consequence and allele frequency |
| `finngen_R14_exome.ld_annotated.tsv.gz` | Exome variants in LD with FinnGen credible set leads, annotated with phenotype and credible set metadata (~900k rows) |
| `plink_fg_merged_chr/finngen_R14_exome_chr[N].{bed,bim,fam}` | Per-chromosome merged plink dataset (FG array + exome variants, 43,289 samples, chr1–23); variant IDs in `CHROM_POS_REF_ALT` format, ~31.1 M variants genome-wide |

### Documentation

#### Logs

| File | Description |
|---|---|
| `[EXOME_DATASET].QC_ANNOTATED.report.txt` | Per-dataset QC filtering statistics by chromosome |
| `hwe/[EXOME_DATASET].hwe_summary.tsv.gz` | Per-dataset Hardy-Weinberg equilibrium + allele frequency summary, autosomes only |

#### Figures

| File | Description |
|---|---|
| `FG_EXOME_resolved_flowchart.png` | Flowchart of sample ID resolution across datasets |
