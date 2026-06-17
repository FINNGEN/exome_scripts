# EXOME data processing

Scripts and WDL workflows for QC-filtering and sample-matching multiple exome cohorts (BOTNIA, ADPKD, DALY, WES) against the FinnGen R14 reference panel. The primary goal is to identify which exome samples correspond to FinnGen participants, enabling downstream data integration. Sample matching uses KING --duplicate on a curated set of ~10,000 HM3 SNPs. Variant QC and genotype filtering are handled by a set of bcftools-based WDL workflows.

---

## Summary results

<!-- BEGIN:data/FG_EXOME_resolved_stats.md -->
## Mapping Totals

| GROUP | TOTAL | ADPKD | BOTNIA | DALY | WES | PCT | NOTES |
| --- | --- | --- | --- | --- | --- | --- | --- |
| MATCHED | 43302 | 619 | 7031 | 12233 | 23419 | 95.4% | samples with a final QRY→REF mapping in the output |
| DROPPED | 1456 | 10 | 22 | 119 | 1305 | 3.2% | found by KING but excluded from final mapping |
| NO MATCH | 636 | 0 | 111 | 53 | 472 | 1.4% | absent from ref or below KING concordance threshold |
| TOTAL | 45394 | 629 | 7164 | 12405 | 25196 | 100.0% |  |

## Mapping Breakdown

| GROUP | STATUS | TOTAL | ADPKD | BOTNIA | DALY | WES | PCT | NOTES |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| MATCHED | ID_CONFIRMED | 40189 | 601 | 5937 | 11816 | 21835 | 88.5% | single candidate; KING match confirmed by matching IDs |
| MATCHED | RESOLVED_BY_ID | 160 | 6 | 31 | 102 | 21 | 0.4% | twins in ref; query ID matched one candidate |
| MATCHED | RESOLVED_BY_ALIAS | 1524 | 0 | 1054 | 0 | 470 | 3.4% | twins in ref; candidates are known aliases of each other |
| MATCHED | UNIQUE | 52 | 0 | 2 | 11 | 39 | 0.1% | single candidate; matched by genetics only |
| MATCHED | CONFLICT_KEPT | 1377 | 12 | 7 | 304 | 1054 | 3.0% | contested ref ID; kept after priority tiebreak; 1377 ref IDs contested, avg 31.4 queries/ref |
| DROPPED | CONFLICT_DROPPED | 1456 | 10 | 22 | 119 | 1305 | 3.2% | contested ref ID; lost tiebreak; REF_MAPPED = NA |
| NO MATCH | MISSING | 636 | 0 | 111 | 53 | 472 | 1.4% | no KING match found |
<!-- END:data/FG_EXOME_resolved_stats.md -->

### Mapping flowchart

![Mapping flowchart](data/FG_EXOME_resolved_flowchart.png)



## SAMPLE MATCHING

Sample matching identifies which exome samples correspond to samples in the FinnGen plink reference panel using KING kinship (`exome_duplicates.wdl`), followed by post-processing with `scripts/resolve_mapping.py` to produce a clean final mapping.

### exome_duplicates.wdl

KING-based duplicate detection between exome VCFs and the FinnGen plink reference. Preferred over gtcheck when sample counts are large — KING scales better and gives a cleaner kinship coefficient rather than a discordance rate.

**How the workflow works:**

```
MakeRegionSnplists
      │
      └─ FilterVCF ×(n_vcfs × n_regions)   [parallel scatter]
             │
       ConcatVCF                            [per VCF]
             │
       SubsetQuery (VCF → plink, shared HM3 variants)
             │
        ┌────┴────────────────┐
   SubsetRef              FilterSNPs
   (ref plink →           (QC filter query plink,
    shared variants)       build final SNP list)
        │                      │
        └────────┬─────────────┘
            PrepQuery / PrepRef
            (subset to QC snplist, het-filter,
             tag IDs, split into chunks)
                   │
              KingShards
              (KING --duplicate, all chunk pairs)
                   │
            SummarizeKing
            (per-sample duplicate summary)
                   │
             GatherResults
             (combine summaries, plots, global stats,
              resolve ambiguities, final mapping + stats)
```

**Step-by-step:**

1. **MakeRegionSnplists**: Downloads Berisa LD blocks (or takes a user-supplied file), assigns HM3 SNPs to blocks, and greedily merges blocks into `n_regions` roughly equal-sized regions. Outputs one SNP list per region (zero-padded filenames so glob order is numeric).

2. **FilterVCF** *(scatter over all VCF × region combinations)*: For each region, runs `bcftools view -R/-T` against the remote GCS VCF using a tabix pre-filter for speed. Outputs one small `vcf.gz` per shard.

3. **ConcatVCF** *(per VCF)*: Sorts the shard VCFs by region index (extracted from filename, not shard order), concatenates with `bcftools concat`, and indexes the result. The `.tbi` is a declared output so it is cached and available for downstream tasks.

4. **SubsetQuery**: Converts the concatenated VCF to plink1 binary (`--make-bed`) extracting only variants present in the HM3 bim. Handles VCF input by detecting the `.vcf.gz` extension — no `input_type` flag needed.

5. **SubsetRef** *(runs in parallel with FilterSNPs)*: Subsets the reference plink to the same shared variants, using the query bim as the SNP list (the task normalises bim → ID column automatically).

6. **FilterSNPs** *(runs in parallel with SubsetRef)*: Applies QC filters to the query plink (`--geno`, `--hwe`, `--maf`). High-quality variants go to the top of the list; the remaining variants are shuffled and appended below. `head -n target_snps | sort -V` then takes the top N sorted by genomic position. Outputs both the final SNP list and the HQ-only list for parameter inspection.

7. **PrepQuery / PrepRef**: Subsets each plink dataset to the final SNP list, removes samples with inbreeding coefficient F > `max_het_F`, annotates IIDs with a dataset prefix (`PREFIX_SAMPLEID`) so KING can distinguish query from ref samples, and splits into chunks of `chunk_size` samples for parallel KING. Note: only IID is prefixed, not FID — the prefix is used internally and stripped back when building the summary.

8. **KingShards**: Runs `king --duplicate` across all query-chunk × ref-chunk pairs sequentially within the task. Filters the `.con` output to only cross-dataset pairs (one sample has the query prefix on IID, the other does not), merges, and gzips the result.

9. **SummarizeKing**: Uses FID (never prefixed) to join the merged `.con.gz` against the query `.fam`, producing a per-sample TSV: each row is one query sample with a comma-separated list of matching reference FIDs, or `MISSING` if none found. Also generates a concordance diagnostic PNG (concordance distribution, IBS0 vs concordance scatter, SNP count per pair).

10. **GatherResults** *(always runs, even if some datasets fail)*: Collects all per-dataset summaries, adds a `DATASET` column (query prefix only), and concatenates into an intermediate `{plink_prefix}_EXOME_summary.tsv`. Stacks concordance PNGs into `{plink_prefix}_EXOME_concordance.png`. Then runs the resolve_mapping logic — with optional alias file for twin disambiguation — to produce `{plink_prefix}_EXOME_resolved.tsv`, `_stats.tsv`, and `_stats.md`.

**Inputs:**

```json
{
  "exome_duplicates.vcf_pairs":    [["ADPKD", "gs://bucket/adpkd.vcf.gz"]],
  "exome_duplicates.plink_bed":    "gs://bucket/finngen_R14_hm3.bed",
  "exome_duplicates.plink_prefix": "FG",
  "exome_duplicates.aliases":      "gs://bucket/finngen_R14_duplicate_list.txt",
  "exome_duplicates.n_regions":    100,
  "exome_duplicates.target_snps":  10000,
  "exome_duplicates.max_het_F":    0.3,
  "exome_duplicates.chunk_size":   10000
}
```

`berisa_blocks` and `aliases` are optional. If `aliases` is omitted, alias-based resolution is skipped. The Berisa LD block file is downloaded automatically if not supplied.

**Outputs:**

- `subset_vcfs[]` / `subset_tbis[]`: Concatenated HM3-subset VCF + index per input dataset (cached for reuse)
- `snplists[]`: Final SNP list used for KING (HQ variants + random padding, sorted by position)
- `hq_snplists[]`: HQ-only SNP list before padding — inspect to tune QC thresholds
- `duplicates_con[]`: Gzipped KING `.con` file with all cross-dataset duplicate pairs
- `summary[]`: Per-sample TSV — one query sample per row, matched reference IDs or `MISSING`
- `concordance_plots[]`: Per-dataset concordance diagnostic PNGs
- `excluded_samples_query[]` / `excluded_samples_ref[]`: Het-outlier samples removed before KING
- `combined_plot`: `{plink_prefix}_EXOME_concordance.png` — all concordance PNGs stacked
- `resolved_mapping`: `{plink_prefix}_EXOME_resolved.tsv` — final QRY→REF mapping with columns `QUERY`, `REF_MAPPED`, `DATASET`, `STATUS`, `CANDIDATES`, `ALIAS_NOTE`
- `resolved_stats_tsv`: `{plink_prefix}_EXOME_resolved_stats.tsv` — group totals + per-status breakdown with per-dataset counts
- `resolved_stats_md`: `{plink_prefix}_EXOME_resolved_stats.md` — same stats in Markdown, ready to paste into this README

---

### Resolve mapping logic

Runs inside `GatherResults` on `combined_summary.tsv`. Also available standalone via `scripts/resolve_mapping.py` for local re-runs without re-executing the full WDL.

```bash
python scripts/resolve_mapping.py combined_summary.tsv [--out my_mapping.tsv] [--seed 42]
```

Handles two real-world complications: twins in the reference (one query matches multiple ref candidates) and true duplicates within a query cohort (multiple queries match the same ref ID).

**Two passes:**

1. **Categorise** — each row classified independently:
   - Single candidate, query ID starts with the ref ID → `ID_CONFIRMED` (genetic + ID agreement, strongest evidence); handles exact matches and suffix variants (e.g. `FGXXXXX_dup1` → ref `FGXXXXX`) without assuming a specific suffix format
   - Single candidate, IDs differ → `UNIQUE` (genetics only, the normal case)
   - Multiple candidates, query ID is one of them → `RESOLVED_BY_ID` (twins in ref; ID identifies the right one)
   - Multiple candidates, query is a known alias of one candidate → `RESOLVED_BY_ALIAS`
   - Multiple candidates, no resolution → `AMBIGUOUS_UNRESOLVED`

2. **Surjectivity check** — each REF ID must appear at most once in the final mapping. Conflicts (multiple queries claiming the same ref) are broken by random draw for now (`CONFLICT_KEPT[orig]` / `CONFLICT_DROPPED[orig]`). TODO: replace with QC tiebreaker (concordance score, het-F, n_snps).

**Output groups:**

| Group | Statuses | Meaning |
|-------|----------|---------|
| **MATCHED** | `ID_CONFIRMED`, `RESOLVED_BY_ID`, `RESOLVED_BY_ALIAS`, `UNIQUE`, `CONFLICT_KEPT[...]` | Has a final QRY→REF mapping in the output |
| **DROPPED** | `CONFLICT_DROPPED[...]`, `AMBIGUOUS_UNRESOLVED` | Found by KING but excluded from final mapping |
| **NO MATCH** | `MISSING` | No KING match found |

**Alias handling** — an optional tab-delimited file lists known alias groups (one group per line, space/tab-separated). For each REF candidate, the script checks whether the query starts with any member of that candidate's alias group. This naturally handles suffix variants (e.g. `FGXXXXX_dup1`) without requiring them to be registered in the alias file. A single matching candidate → `RESOLVED_BY_ALIAS`; multiple matching candidates → `AMBIGUOUS_UNRESOLVED`. Alias IDs cannot appear in `REF_MAPPED`; a post-processing check enforces this.

**Test mode** — a self-contained test dataset covering every status category can be run without any input files:

```bash
python scripts/resolve_mapping.py --test
```

The tables below show the built-in input, alias groups, resolved mapping, and check results (auto-updated by `scripts/resolve_mapping.py --test`):

<!-- BEGIN:test/output.md -->
### Input

| DATASET | QUERY | DUPLICATES |
| --- | --- | --- |
| ds1 | REF001 | REF001 |
| ds1 | REF010_v2 | REF010 |
| ds1 | QRY001 | REF002 |
| ds1 | QRY_ALIAS | REF_ALIAS_A |
| ds1 | QRY_ALIAS_B_v2 | REF_ALIAS_B |
| ds1 | REF_T1 | REF_T1,REF_T2 |
| ds1 | QRY_T_ALIAS | REF_T2,REF_T3 |
| ds2 | QRY003 | REF003,REF004 |
| ds2 | QRY004 | REF003 |
| ds2 | QRY005 | REF005,REF006 |
| ds2 | QRY006 | REF005 |
| ds2 | QRY007 | REF006 |
| ds2 | QRY008 | MISSING |
| ds2 | QRY009 | REF007 |
| ds2 | QRY010 | REF007 |
| ds2 | QRY011 | REF008,REF009 |

### Aliases

| Group |
| --- |
| `QRY_ALIAS` ↔ `REF_ALIAS_A` |
| `QRY_ALIAS_B` ↔ `REF_ALIAS_B` |
| `QRY_T_ALIAS` ↔ `REF_T2` |

### Output

| DATASET | QUERY | DUPLICATES | REF_MAPPED | STATUS | ALIAS_NOTE |
| --- | --- | --- | --- | --- | --- |
| ds1 | REF001 | REF001 | REF001 | ID_CONFIRMED[ID_CONFIRMED] | — |
| ds1 | REF010_v2 | REF010 | REF010 | ID_CONFIRMED[ID_CONFIRMED] | — |
| ds1 | QRY001 | REF002 | REF002 | UNIQUE[UNIQUE] | — |
| ds1 | QRY_ALIAS | REF_ALIAS_A | REF_ALIAS_A | RESOLVED_BY_ALIAS[RESOLVED_BY_ALIAS] | — |
| ds1 | QRY_ALIAS_B_v2 | REF_ALIAS_B | REF_ALIAS_B | RESOLVED_BY_ALIAS[RESOLVED_BY_ALIAS] | — |
| ds1 | REF_T1 | REF_T1,REF_T2 | REF_T1 | RESOLVED_BY_ID[RESOLVED_BY_ID] | — |
| ds1 | QRY_T_ALIAS | REF_T2,REF_T3 | REF_T2 | RESOLVED_BY_ALIAS[RESOLVED_BY_ALIAS] | — |
| ds2 | QRY003 | REF003,REF004 | AMBIGUOUS | AMBIGUOUS_UNRESOLVED[AMBIGUOUS_UNRESOLVED] | query_not_in_alias_file |
| ds2 | QRY004 | REF003 | REF003 | UNIQUE[UNIQUE] | — |
| ds2 | QRY005 | REF005,REF006 | AMBIGUOUS | AMBIGUOUS_UNRESOLVED[AMBIGUOUS_UNRESOLVED] | query_not_in_alias_file |
| ds2 | QRY006 | REF005 | REF005 | UNIQUE[UNIQUE] | — |
| ds2 | QRY007 | REF006 | REF006 | UNIQUE[UNIQUE] | — |
| ds2 | QRY008 | MISSING | NA | MISSING[MISSING] | — |
| ds2 | QRY009 | REF007 | REF007 | CONFLICT_KEPT[UNIQUE] | — |
| ds2 | QRY010 | REF007 | NA | CONFLICT_DROPPED[UNIQUE] | — |
| ds2 | QRY011 | REF008,REF009 | AMBIGUOUS | AMBIGUOUS_UNRESOLVED[AMBIGUOUS_UNRESOLVED] | query_not_in_alias_file |

### Checks

| QUERY | STATUS | REF_MAPPED |  |
| --- | --- | --- | --- |
| REF001 | ID_CONFIRMED[ID_CONFIRMED] | REF001 | ✓ |
| REF010_v2 | ID_CONFIRMED[ID_CONFIRMED] | REF010 | ✓ |
| QRY001 | UNIQUE[UNIQUE] | REF002 | ✓ |
| QRY_ALIAS | RESOLVED_BY_ALIAS[RESOLVED_BY_ALIAS] | REF_ALIAS_A | ✓ |
| QRY_ALIAS_B_v2 | RESOLVED_BY_ALIAS[RESOLVED_BY_ALIAS] | REF_ALIAS_B | ✓ |
| REF_T1 | RESOLVED_BY_ID[RESOLVED_BY_ID] | REF_T1 | ✓ |
| QRY_T_ALIAS | RESOLVED_BY_ALIAS[RESOLVED_BY_ALIAS] | REF_T2 | ✓ |
| QRY003 | AMBIGUOUS_UNRESOLVED[AMBIGUOUS_UNRESOLVED] | AMBIGUOUS | ✓ |
| QRY004 | UNIQUE[UNIQUE] | REF003 | ✓ |
| QRY005 | AMBIGUOUS_UNRESOLVED[AMBIGUOUS_UNRESOLVED] | AMBIGUOUS | ✓ |
| QRY006 | UNIQUE[UNIQUE] | REF005 | ✓ |
| QRY007 | UNIQUE[UNIQUE] | REF006 | ✓ |
| QRY008 | MISSING[MISSING] | NA | ✓ |
| QRY011 | AMBIGUOUS_UNRESOLVED[AMBIGUOUS_UNRESOLVED] | AMBIGUOUS | ✓ |
| QRY009+QRY010 | 1×CONFLICT_KEPT + 1×CONFLICT_DROPPED | — | ✓ |

**15/15 checks — all passed**
<!-- END:test/output.md -->

---

### exome_rename.wdl

**Sample subsetting and renaming using the QRY→REF mapping from `exome_duplicates.wdl`**

**What it does:**

For each dataset × chromosome, streams the QC-annotated VCF from GCS in size-based chunks, subsets to matched samples, and renames sample IDs from query IDs to FinnGen IDs. Outputs one VCF per dataset × chromosome.

```
QueryChromPositions[C]   one task per chrom — queries variant positions
                         for all datasets in parallel (output is cached)
        │
BuildAllRegions          single task — derives size-based chunk regions
                         per dataset × chrom from position files
        │
SubsetChunk[D×C×chunks]  scatter — streams one genomic chunk from GCS,
                         subsets to matched samples, renames to FinnGen IDs
        │
ConcatChromVCF[D×C]      scatter — selects matching chunks and concatenates
                         into one VCF per dataset + chrom
```

**Inputs:**

```json
{
  "exome_rename.resolved_mapping": "gs://bucket/FG_EXOME_resolved.tsv",
  "exome_rename.vcf_pairs":        [["BOTNIA", "gs://bucket/BOTNIA.QC_ANNOTATED.vcf.gz"], ...],
  "exome_rename.chunk_mb":         500,
  "exome_rename.suffix":           "fg_ids"
}
```

`resolved_mapping` is the `_resolved.tsv` output from `exome_duplicates.wdl` — only rows with a non-NA, non-AMBIGUOUS `REF_MAPPED` are used. `vcf_pairs` defaults to the four project datasets if omitted.

**Outputs:**

- `chrom_vcfs[]`: One VCF per dataset × chromosome, samples renamed to FinnGen IDs, filename pattern `{base}.QC_ANNOTATED_{suffix}_{chrom}.vcf.gz`
- `chrom_tbis[]`: Corresponding index files

---

## LD step

The LD step computes linkage disequilibrium (LD) between FinnGen array variants and exome-only variants — i.e. variants present in the exome cohorts but absent from the FinnGen plink reference. This is used to identify which exome variants are well-tagged by FinnGen array variants, enabling downstream imputation and fine-mapping analyses.

### exome_ld.wdl

**Merge FG + exome plink data per chromosome, compute LD, flag coding status, and gather summary**

**Workflow phases:**

```
BuildFgRegions
      │   estimates chunk sizes from FG VCF index + exome sample fraction;
      │   splits per-chrom BIM positions into chunk_mb-sized regions
      │
SubsetFgChunk ×(n_chroms × n_chunks)   [parallel scatter]
      │   streams each FG VCF region from GCS FUSE, subsets to exome union samples
      │
ConcatFgChrom ×n_chroms
      │   concatenates chunks → one VCF.gz per chromosome
      │
VcfToPlink (FG) ×n_chroms             [runs in parallel with Concat]
      │   FG subsetted VCF → plink1 BED, updates sex from phenotype file
      │
VcfToPlink (exome) ×(n_datasets × n_chroms)
      │   each exome VCF → plink1 BED, excluding variants already in FG BIM
      │
MergeChrom ×n_chroms
      │   plink --merge-list: FG bed + all exome beds → one merged BED per chrom
      │
      ├─ ComputeLd
      │     plink2 --r2-unphased anchored on FG variants (--ld-snp-list fg.bim)
      │     outputs native .vcor.zst (no custom docker)
      │
      ├─ FilterLd (calls flag_ld_coding.py)
      │     decompresses .vcor.zst → .vcor.gz; keeps FG→exome pairs with R2≥min_r2;
      │     optionally adds is_fg_coding/is_ex_coding
      │     outputs exome_finngen_ld_<chrom>.ld.tsv.gz + .vcor.gz
      │
      └─ PlinkToVcf
            exports merged plink → bgzipped VCF (used as input to VEP annotation)

GatherLd   [after scatter]
      concatenates all per-chrom ld.tsv.gz and vcor.gz files;
      runs summarize_ld.py → stats TSV + 4 summary figures
```

**Step-by-step:**

1. **BuildFgRegions**: Queries the FG VCF index (no download) to estimate bytes-per-variant and the sample-count ratio between the exome union and the full FG cohort. Uses these to derive a target number of variants per chunk so that each chunk produces approximately `chunk_mb` MB of output. Splits each chromosome's position list from `positions_bim` accordingly and writes a task table with three columns: `fg_vcf_path | chunk_prefix | region`.

2. **SubsetFgChunk** *(scatter over all chunks)*: Streams one genomic region from the FG VCF via GCS FUSE (`bcftools view -r / -t`), keeping only the union of exome samples. Retries on transient GCS errors.

3. **ConcatFgChrom** *(scatter over chromosomes)*: Sorts chunks by their numeric prefix and concatenates them into a single indexed VCF per chromosome.

4. **VcfToPlink (FG)** *(scatter over chromosomes)*: Converts each FG per-chrom VCF to plink1 BED format using plink2. Updates sample sex from the FinnGen phenotype file.

5. **VcfToPlink (exome)** *(scatter over datasets × chromosomes)*: Converts each exome VCF to plink1 BED, with `--exclude` on the FG BIM so the resulting fileset contains only exome-private variants.

6. **MergeChrom** *(scatter over chromosomes)*: Merges the FG plink fileset and all exome plink filesets for a given chromosome into a single combined BED via `plink --merge-list`.

7. **ComputeLd** *(scatter over chromosomes)*: Runs `plink2 --r2-unphased zs` on the merged BED using the FG BIM as `--ld-snp-list`, so LD is computed only for pairs where one variant is a FG array variant. The `zs` modifier tells plink2 to write the native `ld.vcor.zst` (zstd-compressed); `--zst-level 1` sets the compression level. No custom docker — runs on the same VM as the other plink tasks.

8. **FilterLd** *(scatter over chromosomes)*: Decompresses the `.vcor.zst` to `.vcor.gz` (bgzip), then calls `scripts/flag_ld_coding.py` to keep only FG→exome pairs with R² ≥ `min_r2` (default 0.6) and optionally annotate coding status. Outputs both `exome_finngen_ld_<chrom>.ld.tsv.gz` (filtered pairs) and `exome_finngen_ld_<chrom>.vcor.gz` (all FG-anchored pairs above the plink2 window threshold).

9. **PlinkToVcf** *(scatter over chromosomes)*: Exports the merged plink BED back to a bgzipped VCF. These VCFs are the input for the VEP annotation step described below.

10. **GatherLd**: Concatenates all per-chrom `ld.tsv.gz` and `vcor.gz` files (header kept from first chrom, skipped for rest) into genome-wide combined files. Then runs `scripts/summarize_ld.py` on the per-chrom LD files to produce a stats TSV and four summary figures (see below).

**Inputs:**

```json
{
  "exome_ld.fg_vcf_template":  "gs://bucket/finngen_R14_chrCHROM.vcf.gz",
  "exome_ld.fg_pheno_file":    "gs://bucket/finngen_R14_minimum_1.0.txt.gz",
  "exome_ld.exome_vcf_pairs":  [["ADPKD", "gs://bucket/adpkd_chrCHROM.vcf.gz"], ...],
  "exome_ld.positions_bim":    "gs://bucket/finngen_R14.bim",
  "exome_ld.chroms":           ["1","2",...,"22","23"],
  "exome_ld.chunk_mb":         600,
  "exome_ld.ld_params":        "--ld-window-kb 1000 --ld-window-r2 0.05",
  "exome_ld.min_r2":           0.6,
  "exome_ld.out_prefix":       "finngen_R14_exome",
  "exome_ld.annot":            "gs://bucket/vep_annotation.pkl",
  "exome_ld.exome_docker":     "eu.gcr.io/finngen-refinery-dev/exome_bioinf:ld"
}
```

`annot` is optional. `CHROM` in VCF template paths is replaced at runtime with each chromosome name. `min_r2` controls the R² threshold applied in `FilterLd` (default 0.6); note this is separate from `--ld-window-r2` in `ld_params` which is the lower bound passed to plink2 (default 0.05).

**Outputs:**

- `merged_plink[][]`: Per-chromosome merged plink filesets (BED/BIM/FAM/log) — FG + all exome datasets combined
- `ld_results[]`: Per-chrom `exome_finngen_ld_<chrom>.ld.tsv.gz` — FG→exome pairs with R²≥min_r2, with optional coding flags
- `vcor_results[]`: Per-chrom `exome_finngen_ld_<chrom>.vcor.gz` — all FG-anchored pairs above the plink2 window threshold
- `merged_vcf[]`: Per-chromosome bgzipped VCFs of the merged plink data (used as input to VEP)
- `ld_combined`: `{out_prefix}.ld.tsv.gz` — all chromosomes concatenated
- `vcor_combined`: `{out_prefix}.vcor.gz` — all chromosomes concatenated
- `ld_stats`: `{out_prefix}_ld_stats.tsv` — per-chromosome summary table (see `summarize_ld.py`)
- `figures[]`: Four summary PNGs (see `summarize_ld.py`)

---

### flag_ld_coding.py

`scripts/flag_ld_coding.py` is called by the `FilterLd` task. It can also be run standalone for local re-runs.

**What it does:**

1. Loads the FG variant IDs from the FG BIM file.
2. Reads the raw plink2 `.vcor` file (plain or gzipped).
3. Drops pairs with R² below `--min_r2` (default 0.6).
4. Keeps only rows where `ID_A` is a FG variant and `ID_B` is not — i.e. FG→exome pairs. This drops exome→exome and FG→FG pairs that plink2 may include.
5. If `--annot` is provided, looks up the most severe VEP consequence for each variant and adds boolean columns `is_fg_coding` / `is_ex_coding` (true if the consequence is missense, stop-gain, frameshift, splice, or similar protein-altering change).
6. On first read the annotation TSV is parsed and cached as a `.pkl` file in the working directory. Subsequent runs (including the WDL task) load the `.pkl` directly, which is much faster than re-reading the full TSV.

**Output columns:**

| Column | Description |
|--------|-------------|
| `FG_SNP` | FinnGen array variant ID |
| `EXOME_SNP` | Exome-private variant ID |
| `R2` | Unphased r² between the pair |
| `is_fg_coding` | `True` if FG variant has a coding VEP consequence *(only with --annot)* |
| `is_ex_coding` | `True` if exome variant has a coding VEP consequence *(only with --annot)* |

**Annotation file format:**

The annotation file produced by `run_vep_annotate.sh` is a bgzipped TSV (`exome_sites_only_annotated_annot.tsv.bgz`) with the following columns:

```
1  locus
2  alleles
3  rsid
4  variant            ← used as the variant ID key (format: CHROM:POS:REF:ALT)
5  gene_most_severe
6  most_severe        ← VEP consequence used to determine coding status
7  genes_most_severe
```

The script reads `variant` and `most_severe`, builds a `variant → consequence` dictionary, and saves it as a `.pkl`. The `variant` column must match the normalised ID format used in the `.vcor` file (`chr`-prefix stripped, underscores replaced with colons).

**Standalone usage:**

```bash
python scripts/flag_ld_coding.py \
  --vcor    exome_finngen_ld_21.vcor.gz \
  --fg_bim  fg_chr21.bim \
  --annot   exome_sites_only_annotated_annot.tsv.bgz \
  --min_r2  0.6 \
  --out     exome_finngen_ld_21.ld.tsv
```

Run this locally once with the full annotation TSV — this generates the `.pkl` cache file in the working directory. Upload that `.pkl` to GCS and pass it as `exome_ld.annot` in Pass 2 so the WDL task loads it directly without re-parsing the TSV.

---

### summarize_ld.py

`scripts/summarize_ld.py` is called by the `GatherLd` task and can also be run locally on a set of per-chrom LD files.

**What it does:**

Reads N per-chromosome `ld.tsv.gz` files, computes per-chromosome statistics on unique variants and pairs, and writes four summary figures plus a stats TSV.

**Outputs:**

| File | Description |
|------|-------------|
| `{prefix}_ld_stats.tsv` | Per-chromosome counts: n_pairs, n_fg/ex_variants, coding fractions, pair breakdown |
| `{prefix}_fig1_variants.png` | Stacked bar: unique coding/non-coding variants per chrom (FG and exome panels) |
| `{prefix}_fig2_pairs.png` | Stacked bar: pairs by coding category per chrom (both/FG-only/exome-only/neither) |
| `{prefix}_fig3_r2_dist.png` | Violin: R² distribution by coding category, pooled across all chroms |
| `{prefix}_fig4_coding_frac.png` | Line plot: coding fraction per chrom for FG and exome variants |

Stats TSV columns: `chrom`, `n_pairs`, `n_fg_variants`, `fg_coding`, `fg_coding_pct`, `n_ex_variants`, `ex_coding`, `ex_coding_pct`, `n_pairs_both_coding`, `n_pairs_fg_only`, `n_pairs_ex_only`, `n_pairs_neither`.

**Standalone usage:**

```bash
# from a directory containing the per-chrom files
python scripts/summarize_ld.py \
  <(ls exome_finngen_ld_*.ld.tsv.gz) \
  --prefix finngen_R14_exome \
  --outdir ~/results/ld_summary

# or pipe from stdin
ls exome_finngen_ld_*.ld.tsv.gz | python scripts/summarize_ld.py - --outdir .

# quick test on first 1000 lines per file
python scripts/summarize_ld.py <(ls *.ld.tsv.gz) --test
```

---

### VEP annotation

To add coding status flags to the LD output, VEP annotations must be generated from the merged VCFs produced by `exome_ld.wdl`. The recommended approach is a two-pass run:

**Pass 1 — run without annotation:**

Run `exome_ld.wdl` without the `annot` input. The workflow produces `merged_vcf[]` outputs (one per chromosome) alongside the unannotated `ld_results[]`.

**Annotate the VCFs:**

Run `run_vep_annotate.sh` from [FinnGen commons](https://github.com/FINNGEN/commons/blob/master/variant_annotation/scripts/run_vep_annotate.sh). The script does not take VCF paths as arguments — instead, edit the `vcf_in` variable near the top of the script to point to the `merged_vcf[]` outputs from Pass 1:

```bash
# inside run_vep_annotate.sh, change this line:
vcf_in="gs://r12-data/exome/renamed_final/vcf/*.gz"
# to the GCS path of the VCFs produced by PlinkToVcf, e.g.:
vcf_in="gs://your-bucket/finngen_R14_exome_chr*.vcf.gz"
```

Then run the script:

```bash
bash commons/variant_annotation/scripts/run_vep_annotate.sh
```

This produces a VEP annotation TSV (or pickle) containing variant IDs and `most_severe` consequence columns consumed by `flag_ld_coding.py`.

**Pass 2 — re-run with annotation:**

Before re-submitting the WDL, run `flag_ld_coding.py` locally once on any single-chromosome `.vcor` file with the annotation TSV. This generates the `.pkl` cache (see the standalone usage in the `flag_ld_coding.py` section above). Upload the `.pkl` to GCS, then re-submit `exome_ld.wdl` pointing `annot` at it. The `FilterLd` task loads the `.pkl` directly and the final `ld_results[]` files will include the `is_fg_coding` / `is_ex_coding` columns.

```json
{
  "exome_ld.annot": "gs://bucket/exome_sites_only_annotated_annot.pkl"
}
```

---

## Annotation/QC

This section summarizes the datasets and processing steps required to merge them into a single quality-controlled dataset.

### Datasets

#### Daly Data
- **Samples**: 12,405 samples with non-FG IDs
- **Status**: Merged into a single file by Lea
- **Issues**: Required re-heading due to header formatting problems
- **QC Approach**: Contains FILTER entries allowing quick filtering with `FILTER~"NO_HQ_GENOTYPES"`
- **Notes**: [GATK GVS filter explanations](https://github.com/broadinstitute/gatk/blob/ah_var_store/scripts/variantstore/beta_docs/gvs-outputs.md#applying-the-gvs-joint-calling-filters)

#### ADPKD
- **Samples**: 629 samples with FG IDs
- **Status**: Merged by Lea and released to red sandbox
- **Issues**: Required re-heading (bcftools warnings): `[W::bcf_hdr_check_sanity] PL should be declared as Number=G`'. 
- **QC Approach**: Contains AC, DP, and GQ fields for standard filtering
- **Sample Matching**:

| Match Type | Count |
|------------|-------|
| Direct matches | 629 |
| Indirect matches via mapping | 0 |
| No match | 1 |

#### BOTNIA
- **Samples**: 7,164 samples with FG IDs
- **Status**: Single file, ready for processing
- **QC Approach**: Contains AC, DP, and GQ fields for standard filtering
- **Sample Matching**:

| Match Type | Count |
|------------|-------|
| Direct matches | 5,973 |
| Indirect matches via mapping | 1,081 |
| No match | 110 |


#### FINNGEN-WES
25201 with FGID (Sometimes with "_dup") extension to indicate duplicate. In that case there are only 23223 unique IDs left

- **Samples**:25201 with FGID (Sometimes with "_dup") extension to indicate duplicate. In that case there are only 23223 left
- **Status**: Chrom files
- **Issues** Chrom files are *massive* and mostly empty.
- **QC Approach**: Contains AC, DP, and GQ fields for standard filtering
- **Sample Matching**:

| Match Type | Count |
|------------|-------|
| Direct matches | 22163 unique |
| Indirect matches via mapping | 730 |
| No match | 330 |


### WDL Workflows

This repository contains WDL (Workflow Description Language) workflows for processing VCF files from sequencing data. All workflows are located in the `wdl/` directory.

#### Overview

| Workflow | Purpose | Input Data Type |
|----------|---------|-----------------|
| `wes_chrom.wdl` | **Multi-chromosome parallel filtering** | Per-chromosome VCFs (exome or targeted sequencing) |
| `single_file_qc.wdl` | **Whole genome filtering** | Single whole-genome VCF files (splits by chromosome internally) |
| `daly_qc.wdl` | **Simple filter-based QC** | Any VCF with FILTER flags to remove |

---

#### wes_chrom.wdl

**Parallel filtering workflow for per-chromosome WES data with position-based chunking**

**What it does:**

1. Optionally subsets samples for testing (if `test_sample_count` is provided)
2. Computes statistics on original VCFs (variant counts per chromosome)
3. Pre-filters each VCF (`PreFilter` task):
   - Removes AC=0 variants (monomorphic sites)
   - Annotates variant IDs as `CHROM_POS_REF_ALT`
4. Parallel filters by region using position-based chunking:
   - Splits each chromosome into equal chunks by variant positions
   - **Splits multiallelics** (`bcftools norm -m -any`) and **normalises against reference FASTA** (`bcftools norm -c x`) — done first before any other filters
   - Applies genotype filters (sets low-quality genotypes to missing)
   - Recalculates AC (allele count) after genotype filtering
   - Applies variant filters (removes variants with AC=0, etc.)
5. Validates filtering on sample VCFs (checks filters worked correctly)
6. Creates summary statistics showing variant counts and drop rates per chromosome

**Key features:**

- **Uses both `-r` and `-T` flags**: Fast index-based seeking (`-r`) combined with exact position filtering (`-T`) to avoid overlapping variants between chunks
- **No sorting needed**: Position-based chunking ensures chunks don't overlap, so naive concatenation produces sorted output
- **Extensive validation**: Creates sample VCFs and validates that filters were applied correctly

**Inputs:**

```json
{
  "vcf_list": "path/to/vcf_list.txt",
  "genotype_filter": "FORMAT/DP<10 | FORMAT/GQ<20",
  "variant_filter": "AC>0 & ALT!=\"*\"",
  "cpu_count": 16,
  "norm_fasta": "gs://bucket/reference.fa",
  "test_sample_count": 10  // Optional: subset to first N samples for testing
}
```

Input file format (`vcf_list.txt`):
```
gs://bucket/chr1.vcf.gz
gs://bucket/chr2.vcf.gz
gs://bucket/chr3.vcf.gz
...
```

**Outputs:**

- `filtered_vcfs[]`: Per-chromosome filtered VCF files
- `filtered_vcf_tbis[]`: Corresponding index files
- `merged_vcf`: All chromosomes merged into single VCF
- `merged_vcf_tbi`: Index for merged VCF
- `summary_table`: TSV with chromosome-level and total drop rates
- `original_stats[]`: Variant counts before filtering
- `filtered_stats[]`: Variant counts after filtering
- `validation_reports[]`: Per-chromosome validation reports


---

#### single_file_qc.wdl

**Per-chromosome parallel filtering for whole genome VCF files**

**What it does:**

Designed for **whole genome VCF files** where all chromosomes are in a single file. The workflow automatically splits processing by chromosome for parallelization.

1. Optionally subsets samples for testing
2. Computes chromosome counts for original VCF
3. Filters each chromosome in parallel:
   - **Auto-detects chromosome naming** — if the VCF uses non-prefixed names (`1`, `22`) instead of `chr1`, `chr22`, chromosomes are permanently renamed to `chr`-prefix in the output so both VCF types produce consistently named output
   - Strips null bytes (`tr -d '\0'`) to handle corrupt FORMAT fields
   - **Splits multiallelics** (`bcftools norm -m -any`) and **normalises against reference FASTA** (`bcftools norm -c x`)
   - Sets low-quality genotypes to missing
   - Recalculates AC
   - Filters variants by expression
   - Annotates variant IDs as `CHROM_POS_REF_ALT`
4. Validates filtering worked correctly
5. Outputs filtered VCF files

**Key differences from wes_chrom.wdl:**

- Processes **whole genome files** (not pre-split by chromosome)
- Uses **chromosome-based parallelization** instead of position chunking
- Works with both **chr-prefixed** (e.g. ADPKD) and **non-prefixed** (e.g. BOTNIA) VCFs — output is always chr-prefixed
- No merging step (outputs remain per-chromosome)
- Simpler workflow for whole genome data

**Inputs:**

```json
{
  "vcf_files": ["gs://bucket/whole_genome.vcf.gz"],
  "genotype_filter": "FORMAT/DP<10 | FORMAT/GQ<20",
  "variant_filter": "AC>0 & ALT!=\"*\"",
  "cpu_count": 8,
  "norm_fasta": "gs://bucket/reference.fa",
  "test_sample_count": 10  // Optional
}
```

**Outputs:**

- `filtered_vcfs[]`: Filtered VCF per input file
- `filtered_vcf_tbis[]`: Index files
- `original_chrom_counts[]`: Variant counts per chromosome (before)
- `filtered_chrom_counts[]`: Variant counts per chromosome (after)
- `validation_reports[]`: Validation reports

**How to run:**

```bash
cat > inputs.json << EOF
{
  "single_file_qc.vcf_files": ["gs://bucket/whole_genome.vcf.gz"],
  "single_file_qc.genotype_filter": "FORMAT/DP<10 | FORMAT/GQ<20",
  "single_file_qc.variant_filter": "AC>0 & ALT!=\\\"*\\\"",
  "single_file_qc.cpu_count": 8
}
EOF

java -jar cromwell.jar run wdl/single_file_qc.wdl -i inputs.json
```

---

#### daly_qc.wdl

**FILTER-based QC with header annotation and sample renaming**

**What it does:**

For each input VCF (scatter over vcf_list):

1. **ComputeStats**: Reads the remote index only — gets chrom, variant count, sample count without downloading the VCF.
2. **AnnotateAndRename**: Fetches the VCF header, injects any missing FILTER definitions (`NO_HQ_GENOTYPES`, `ExcessHet`, `LowQual`, `EXCESS_ALLELES`, `OUTSIDE_OF_TARGETS`), builds a `SAMPLE_ID → FINNGENID_finngen` rename map from the `rename_file`, resolves duplicate new IDs with `_dup`/`_dup2` suffixes, and applies header + rename in a single `bcftools reheader` pass (fast — body is copied verbatim).
3. **ParallelFilter**: Splits the chromosome into `cpu_count × chunk_multiplier` equal position windows, filters each chunk in parallel with GNU parallel (`bcftools view -e filter_expression`), annotates variant IDs as `CHROM_POS_REF_ALT`, and concatenates.
4. **ComputeStats** (again): Variant counts on the filtered VCF.
5. **ValidateFiltering**: Checks no variants matching the filter expression remain, verifies `CHROM_POS_REF_ALT` ID format, and reports rename counts.

Then globally:

6. **SummaryStats**: Per-chromosome variant counts + drop rates, derives output root name from input filenames.
7. **SortAndMerge**: Concatenates per-chromosome filtered VCFs (in vcf_list order) into `{root_name}.QC_ANNOTATED.vcf.gz`.

**Inputs:**

```json
{
  "daly_qc.vcf_list":         "path/to/vcf_list.txt",
  "daly_qc.rename_file":      "path/to/rename.tsv",
  "daly_qc.filter_expression": "FILTER~'NO_HQ_GENOTYPES'",
  "daly_qc.cpu_count":        8,
  "daly_qc.vcf_max_gb":       25,
  "daly_qc.chunk_multiplier": 3
}
```

`rename_file` is a TSV with columns: `FINNGENID_finngen(1)`, `FINNGENID_biobank(2)`, `SAMPLE_ID(3)` — samples are renamed from col 3 to col 1. `filter_expression`, `vcf_max_gb`, and `chunk_multiplier` are optional.

**Outputs:**

- `original_stats[]`: Variant/sample counts before filtering (index-only, fast)
- `filtered_stats[]`: Variant counts after filtering
- `validation_reports[]`: Per-VCF filter check + rename summary
- `merged_vcf`: `{root_name}.QC_ANNOTATED.vcf.gz` — all chromosomes merged, samples renamed to FinnGen IDs
- `merged_vcf_tbi`: Index for merged VCF
- `report`: `{root_name}.QC_ANNOTATED.report.txt` — per-chromosome and total drop rates

**How to run:**

```bash
cat > inputs.json << EOF
{
  "daly_qc.vcf_list":         "vcf_files.txt",
  "daly_qc.rename_file":      "rename.tsv",
  "daly_qc.filter_expression": "FILTER~'NO_HQ_GENOTYPES'",
  "daly_qc.cpu_count":        8
}
EOF

java -jar cromwell.jar run wdl/daly_qc.wdl -i inputs.json
```

---

#### Common Patterns

**Filter Expressions**

*Genotype filters* (applied with `bcftools +setGT`):
- `FORMAT/DP<10 | FORMAT/GQ<20` - Set genotypes to missing if DP<10 OR GQ<20
- `FORMAT/DP<10 & FORMAT/GQ<20` - Set genotypes to missing if DP<10 AND GQ<20
- `FORMAT/DP<10` - Only depth filter
- `FORMAT/GQ<20` - Only quality filter

*Variant filters* (applied with `bcftools view -i`):
- `AC>0` - Keep variants with at least 1 alternate allele
- `ALT!="*"` - Remove spanning deletions
- `AC>0 & ALT!="*"` - Both conditions
- `AC>=2` - Keep variants with at least 2 alternate alleles (removes singletons)

**Testing**

All workflows support `test_sample_count` to subset samples for faster testing:

```json
{
  "test_sample_count": 10  // Only process first 10 samples
}
```

**CPU Count**

Controls parallelization level:
- `wes_chrom.wdl`: Number of chunks per chromosome
- `single_file_qc.wdl`: Number of chromosomes processed in parallel
- `daly_qc.wdl`: Number of genomic regions processed in parallel

Recommended: 8-16 for most use cases.

---

#### Troubleshooting

**"Unsorted positions" error during indexing**

*Cause:* Overlapping variants between chunks (old behavior with `-R` only)

*Solution:* Use `wes_chrom.wdl` which uses both `-r` and `-T` to ensure no overlaps

**Filter expression errors**

*Symptom:* `Error occurred while processing the filter`

*Solution:* Ensure special characters are properly escaped:
```json
"variant_filter": "AC>0 & ALT!=\"*\""  // Correct
"variant_filter": "AC>0 & ALT!=*"      // Wrong - missing quotes
```

**Out of memory errors**

*Solution:* Increase memory in runtime section or reduce `cpu_count` (fewer parallel chunks = less memory)

---

#### Requirements

- bcftools (with +setGT and +fill-tags plugins)
- tabix
- parallel (GNU parallel)
- Python 3 with numpy (for position chunking in wes_chrom.wdl)
- WDL runtime (Cromwell, miniwdl, etc.)
