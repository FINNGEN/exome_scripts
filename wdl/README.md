# WDL Workflows

This directory contains WDL (Workflow Description Language) workflows for processing VCF files from whole exome sequencing data.

## Overview

| Workflow | Purpose | Best For |
|----------|---------|----------|
| `wes_chrom.wdl` | **Multi-chromosome parallel filtering** | Processing multiple chromosome VCFs with genotype and variant filters |
| `single_file_qc.wdl` | **Single-file chromosome filtering** | Processing single VCF files with per-chromosome parallel filtering |
| `daly_qc.wdl` | **Simple filter-based QC** | Removing variants with specific FILTER flags |

---

## wes_chrom.wdl

**Parallel filtering workflow for multi-chromosome WES data with position-based chunking**

### What it does

1. **Optionally subsets samples** for testing (if `test_sample_count` is provided)
2. **Computes statistics** on original VCFs (variant counts per chromosome)
3. **Parallel filters by region** using position-based chunking:
   - Splits each chromosome into equal chunks by variant positions
   - Applies genotype filters (sets low-quality genotypes to missing)
   - Recalculates AC (allele count) after genotype filtering
   - Applies variant filters (removes variants with AC=0, etc.)
   - Annotates variant IDs as `CHROM_POS_REF_ALT`
4. **Validates filtering** on sample VCFs (checks filters worked correctly)
5. **Merges filtered VCFs** across chromosomes
6. **Creates summary statistics** showing variant counts and drop rates per chromosome

### Key features

- **Uses both `-r` and `-T` flags**: Fast index-based seeking (`-r`) combined with exact position filtering (`-T`) to avoid overlapping variants between chunks
- **No sorting needed**: Position-based chunking ensures chunks don't overlap, so naive concatenation produces sorted output
- **Extensive validation**: Creates sample VCFs and validates that filters were applied correctly

### Inputs

```json
{
  "vcf_list": "path/to/vcf_list.txt",
  "genotype_filter": "FORMAT/DP<10 | FORMAT/GQ<20",
  "variant_filter": "AC>0 & ALT!=\"*\"",
  "cpu_count": 16,
  "test_sample_count": 10  // Optional: subset to first N samples for testing
}
```

**Input file format** (`vcf_list.txt`):
```
gs://bucket/chr1.vcf.gz
gs://bucket/chr2.vcf.gz
gs://bucket/chr3.vcf.gz
...
```

### Outputs

- `filtered_vcfs[]`: Per-chromosome filtered VCF files
- `filtered_vcf_tbis[]`: Corresponding index files
- `merged_vcf`: All chromosomes merged into single VCF
- `merged_vcf_tbi`: Index for merged VCF
- `summary_table`: TSV with chromosome-level and total drop rates
- `original_stats[]`: Variant counts before filtering
- `filtered_stats[]`: Variant counts after filtering
- `validation_reports[]`: Per-chromosome validation reports

### How to run

```bash
# Create input JSON
cat > inputs.json << EOF
{
  "wes_chrom.vcf_list": "vcf_files.txt",
  "wes_chrom.genotype_filter": "FORMAT/DP<10 | FORMAT/GQ<20",
  "wes_chrom.variant_filter": "AC>0 & ALT!=\\\"*\\\"",
  "wes_chrom.cpu_count": 16
}
EOF

# Submit to Cromwell
java -jar cromwell.jar run wes_chrom.wdl -i inputs.json

# Or with test mode (10 samples only)
cat > test_inputs.json << EOF
{
  "wes_chrom.vcf_list": "vcf_files.txt",
  "wes_chrom.genotype_filter": "FORMAT/DP<10 | FORMAT/GQ<20",
  "wes_chrom.variant_filter": "AC>0 & ALT!=\\\"*\\\"",
  "wes_chrom.cpu_count": 16,
  "wes_chrom.test_sample_count": 10
}
EOF

java -jar cromwell.jar run wes_chrom.wdl -i test_inputs.json
```

### Tasks breakdown

1. **SubsetSamples** (optional): Subsets VCF to first N samples using column cutting
2. **ComputeStats**: 
   - Extracts all variant positions
   - Counts total variants
   - Creates sample VCF (100 variants) for validation
3. **ParallelFilterByRegion**:
   - Splits positions into equal chunks
   - Creates position files (`CHROM POS` format) and region files (`CHROM:START-END`)
   - Runs filtering in parallel using both `-r` and `-T` for efficiency
   - Concatenates results (no sorting needed)
4. **ValidateFiltering**: 
   - Checks genotypes were set to missing correctly
   - Verifies AC field was recalculated
   - Validates variant ID format
5. **SortAndMerge**: Concatenates per-chromosome VCFs into single merged file
6. **SummaryStats**: Creates summary table with drop rates

---

## single_file_qc.wdl

**Per-chromosome parallel filtering for single VCF files**

### What it does

1. Optionally subsets samples for testing
2. Computes chromosome counts for original VCF
3. Filters each chromosome in parallel:
   - Sets low-quality genotypes to missing
   - Recalculates AC
   - Filters variants by expression
   - Annotates variant IDs
4. Validates filtering worked correctly
5. Outputs filtered VCF files

### Key differences from wes_chrom.wdl

- Processes **each chromosome separately** (not position chunks)
- Uses **simple chromosome-based parallelization** instead of position chunking
- No merging step (outputs remain per-chromosome)
- Simpler workflow for single-file processing

### Inputs

```json
{
  "vcf_files": ["gs://bucket/sample.vcf.gz"],
  "genotype_filter": "FORMAT/DP<10 | FORMAT/GQ<20",
  "variant_filter": "AC>0 & ALT!=\"*\"",
  "cpu_count": 8,
  "test_sample_count": 10  // Optional
}
```

### Outputs

- `filtered_vcfs[]`: Filtered VCF per input file
- `filtered_vcf_tbis[]`: Index files
- `original_chrom_counts[]`: Variant counts per chromosome (before)
- `filtered_chrom_counts[]`: Variant counts per chromosome (after)
- `validation_reports[]`: Validation reports

### How to run

```bash
cat > inputs.json << EOF
{
  "single_file_qc.vcf_files": ["gs://bucket/sample.vcf.gz"],
  "single_file_qc.genotype_filter": "FORMAT/DP<10 | FORMAT/GQ<20",
  "single_file_qc.variant_filter": "AC>0 & ALT!=\\\"*\\\"",
  "single_file_qc.cpu_count": 8
}
EOF

java -jar cromwell.jar run single_file_qc.wdl -i inputs.json
```

---

## daly_qc.wdl

**Simple FILTER-based QC workflow**

### What it does

1. Annotates VCF headers with FILTER definitions
2. Filters variants in parallel by genomic regions
3. Removes variants matching a filter expression (default: `FILTER~"NO_HQ_GENOTYPES"`)
4. Merges filtered chromosomes
5. Annotates variant IDs

### Use case

Designed for filtering variants that have already been flagged in the FILTER column. Simpler than the other workflows - just removes flagged variants without genotype-level filtering.

### Inputs

```json
{
  "vcf_list": "path/to/vcf_list.txt",
  "filter_expression": "FILTER~\"NO_HQ_GENOTYPES\"",
  "cpu_count": 8
}
```

### Outputs

- `filtered_vcfs[]`: Per-chromosome filtered VCFs
- `filtered_vcf_tbis[]`: Index files
- `merged_vcf`: Merged output
- `merged_vcf_tbi`: Merged index
- `total_variant_count`: Total variants in merged file

### How to run

```bash
cat > inputs.json << EOF
{
  "daly_qc.vcf_list": "vcf_files.txt",
  "daly_qc.filter_expression": "FILTER~\\\"NO_HQ_GENOTYPES\\\"",
  "daly_qc.cpu_count": 8
}
EOF

java -jar cromwell.jar run daly_qc.wdl -i inputs.json
```

---

## Common Patterns

### Filter Expressions

**Genotype filters** (applied with `bcftools +setGT`):
- `FORMAT/DP<10 | FORMAT/GQ<20` - Set genotypes to missing if DP<10 OR GQ<20
- `FORMAT/DP<10 & FORMAT/GQ<20` - Set genotypes to missing if DP<10 AND GQ<20
- `FORMAT/DP<10` - Only depth filter
- `FORMAT/GQ<20` - Only quality filter

**Variant filters** (applied with `bcftools view -i`):
- `AC>0` - Keep variants with at least 1 alternate allele
- `ALT!="*"` - Remove spanning deletions
- `AC>0 & ALT!="*"` - Both conditions
- `AC>=2` - Keep variants with at least 2 alternate alleles (removes singletons)

### Testing

All workflows support `test_sample_count` to subset samples for faster testing:

```json
{
  "test_sample_count": 10  // Only process first 10 samples
}
```

### CPU Count

Controls parallelization level:
- `wes_chrom.wdl`: Number of chunks per chromosome
- `single_file_qc.wdl`: Number of chromosomes processed in parallel
- `daly_qc.wdl`: Number of genomic regions processed in parallel

Recommended: 8-16 for most use cases.

---

## Troubleshooting

### "Unsorted positions" error during indexing

**Cause**: Overlapping variants between chunks (old behavior with `-R` only)

**Solution**: Use `wes_chrom.wdl` which uses both `-r` and `-T` to ensure no overlaps

### Filter expression errors

**Symptom**: `Error occurred while processing the filter`

**Solution**: Ensure special characters are properly escaped:
```json
"variant_filter": "AC>0 & ALT!=\"*\""  // Correct
"variant_filter": "AC>0 & ALT!=*"      // Wrong - missing quotes
```

### Out of memory errors

**Solution**: Increase memory in runtime section or reduce `cpu_count` (fewer parallel chunks = less memory)

---

## Version History

- **2026-05-08**: Updated `wes_chrom.wdl` to use position-based filtering with `-r` and `-T` flags
- **Initial**: Created workflows for WES filtering

---

## Requirements

- bcftools (with +setGT and +fill-tags plugins)
- tabix
- parallel (GNU parallel)
- Python 3 with numpy (for position chunking in wes_chrom.wdl)
- WDL runtime (Cromwell, miniwdl, etc.)
