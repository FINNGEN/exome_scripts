# EXOME data processing

## SAMPLE MATCHING

Sample matching identifies which exome samples correspond to samples in the FinnGen plink reference panel. This is done using `exome_duplicates.wdl`, which measures genotype concordance between each exome VCF and the reference using `bcftools gtcheck`.

### exome_duplicates.wdl

**What it does:**

For each exome dataset:

1. **SubsetVCF**: Pre-filters the exome VCF to only the positions present in the reference `.bim` file. Runs one `bcftools view -T` job per chromosome in parallel using the tabix index for fast positional lookup. Chunks are concatenated into a single filtered VCF, sample IDs are reheadered with the dataset prefix (e.g. `BOTNIA_SAMPLE1`) so query samples are unambiguously labelled in the gtcheck output, and a SNP ID list is written for the next step.

2. **PlinkSubset**: Subsets the plink reference to the SNPs from SubsetVCF (`--extract`) and renames all sample IIDs in-place to `plink_prefix_SAMPLE` (e.g. `FG_XX00000001`) so reference samples are unambiguously labelled in the gtcheck output.

3. **RunGtcheck**: Splits the renamed plink fam into chunks of `chunk_size` samples, converts each chunk to a bgzipped VCF with plink2 in parallel (32 CPUs, memory divided equally across jobs), indexes chunk_00 and copies the `.tbi` to all other chunks, then runs `bcftools gtcheck --no-HWE-prob` for each chunk in parallel. Results are merged and summarised: for each query sample the best match, second-best match, average concordance across all others, and the ratio best/average are reported.

**Inputs:**

```json
{
  "exome_duplicates.vcf_pairs": [
    ["BOTNIA", "gs://bucket/botnia.vcf.gz"],
    ["ADPKD",  "gs://bucket/adpkd.vcf.gz"]
  ],
  "exome_duplicates.plink_bed":    "gs://bucket/finngen_R14_hm3.bed",
  "exome_duplicates.plink_prefix": "FG",
  "exome_duplicates.bim":          "gs://bucket/custom.bim",
  "exome_duplicates.chunk_size":   100
}
```

`vcf_pairs` is an array of `[prefix, vcf_path]` pairs. `bim` is optional — if omitted the `.bim` derived from `plink_bed` is used. `chunk_size` controls how many reference samples are processed per parallel plink2/gtcheck job (default 100).

**Outputs:**

- `gtcheck_raw[]`: Raw `bcftools gtcheck` output per dataset
- `gtcheck_summary[]`: TSV summary per dataset with columns `QUERY`, `BEST_MATCH`, `BEST_RATE`, `2ND_MATCH`, `2ND_RATE`, `AVG_OTHERS`, `RATIO`

**Reading the output:**

Each row in the summary TSV is one query (exome) sample. `BEST_MATCH` is the most concordant reference sample. `RATIO` = `BEST_RATE / AVG_OTHERS` — a low ratio (< 0.1) indicates a likely true duplicate. `BEST_RATE` is the per-site discordance rate so lower = more similar.

**Local testing:**

`scripts/test_gtcheck.sh` replicates the RunGtcheck step locally:

```bash
scripts/test_gtcheck.sh <query.vcf.gz> <plink_prefix> [output_prefix] [--parallel N]
```


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

**Simple FILTER-based QC workflow**

**What it does:**

1. Computes statistics on original VCFs (variant counts per chromosome)
2. Annotates VCF headers with FILTER definitions (if missing)
3. Filters variants in parallel by genomic regions
4. Removes variants matching a filter expression (default: `FILTER~"NO_HQ_GENOTYPES"`)
5. Computes statistics on filtered VCFs
6. Merges filtered chromosomes
7. Annotates variant IDs
8. Creates summary statistics showing variant counts and drop rates per chromosome

**Use case:**

Designed for filtering variants that have already been flagged in the FILTER column. Simpler than the other workflows - just removes flagged variants without genotype-level filtering.

**Inputs:**

```json
{
  "vcf_list": "path/to/vcf_list.txt",
  "filter_expression": "FILTER~\"NO_HQ_GENOTYPES\"",
  "cpu_count": 8
}
```

**Outputs:**

- `filtered_vcfs[]`: Per-chromosome filtered VCFs
- `filtered_vcf_tbis[]`: Index files
- `original_stats[]`: Variant counts before filtering
- `filtered_stats[]`: Variant counts after filtering
- `merged_vcf`: Merged output
- `merged_vcf_tbi`: Merged index
- `summary_table`: TSV with chromosome-level and total drop rates

**How to run:**

```bash
cat > inputs.json << EOF
{
  "daly_qc.vcf_list": "vcf_files.txt",
  "daly_qc.filter_expression": "FILTER~\\\"NO_HQ_GENOTYPES\\\"",
  "daly_qc.cpu_count": 8
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

