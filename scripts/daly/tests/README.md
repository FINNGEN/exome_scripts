# Test Scripts

This directory contains test and benchmark scripts for VCF processing workflows.

## Scripts

### compare_methods.sh
Benchmarks two methods of parallel VCF processing:
- **Native bcftools**: Uses bcftools with `--threads` flag
- **Region-split method**: Splits VCF into genomic regions and processes chunks in parallel with GNU parallel

Compares execution time and validates that both methods produce identical output. Takes first N variants from input VCF and times both approaches.

**Usage:** `./compare_methods.sh input.vcf.gz N`

---

### process_vcf_parallel.sh
Parallel VCF processing script that splits a VCF by genomic regions and processes chunks concurrently. Uses Python numpy to create equal-sized genomic intervals, processes each chunk with GNU parallel, then concatenates results.

Applies:
- Filtering to remove variants with `NO_HQ_GENOTYPES` filter
- ID annotation using `%CHROM_%POS_%REF_%ALT` pattern

**Usage:** Hardcoded input file - edit script to change input VCF

---

### reheader.sh
Adds missing FILTER header definitions to VCF files. Recreated from WDL task `AnnotateHeaders`. Checks for missing FILTER headers (NO_HQ_GENOTYPES, ExcessHet, LowQual, EXCESS_ALLELES, OUTSIDE_OF_TARGETS) and adds their descriptions to the VCF header.

Supports test mode (processes only first 1,000 variants) and multi-threaded compression. Optimized to copy original index when only headers are modified.

**Usage:** `./reheader.sh <input_vcf> [test_mode] [threads]`

---

### run_parallel_filter_local.sh
Local test script that replicates the `ParallelFilter` WDL task. Splits VCF into equal genomic regions using binary search to find first/last variant positions, processes chunks in parallel with bcftools filtering, and concatenates results.

Applies custom filter expression (default: removes `NO_HQ_GENOTYPES`) and variant ID annotation.

**Usage:** `./run_parallel_filter_local.sh <input.vcf.gz> [filter_expression] [cpu_count]`

---

### test_header_annotation.sh
One-liner test for fast header annotation. Generates FILTER header definitions and injects them into VCF without using bcftools annotate. Processes only the first 100 variants for quick validation.

Tests an alternative approach to header modification and validates the output is parseable.

**Usage:** `./test_header_annotation.sh input.vcf.gz`

---

### test_index_stats.sh
Validates that `bcftools index -s` provides accurate position information. Compares three methods of finding first/last variant positions:
1. Scanning entire VCF with zcat/awk (slow)
2. Using bcftools index stats (fast)
3. Actual last variant position verification

Benchmarks each method and confirms index stats can be trusted for range determination.

**Usage:** `./test_index_stats.sh input.vcf.gz`

---

### test_positions.sh
Binary search-based script to find actual first and last variant positions for each contig in a VCF. Uses bcftools index to get contig information, then performs efficient binary search to locate the last variant position without scanning the entire file.

Outputs a table showing contig, length, variant count, first position, and last position.

**Usage:** `./test_positions.sh <vcf.gz>`
