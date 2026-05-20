#!/bin/bash

# Recreated from WDL task AnnotateHeaders - identical to WDL version
# Usage: ./reheader.sh <input_vcf> [test_mode] [threads]
#   input_vcf: Path to input VCF file (.vcf.gz)
#   test_mode: Optional. Set to "true" to process only first 10k variants (default: false)
#   threads: Optional. Number of threads (default: auto-detect)

if [ $# -lt 1 ]; then
    echo "Usage: $0 <input_vcf> [test_mode] [threads]"
    echo "  input_vcf: Path to input VCF file (.vcf.gz)"
    echo "  test_mode: Optional. Set to 'true' to process only first 10k variants (default: false)"
    echo "  threads: Optional. Number of threads (default: auto-detect)"
    exit 1
fi

input_vcf="$1"
test_mode="${2:-false}"
threads="${3:-}"

# Auto-detect threads if not specified
if [ -z "$threads" ]; then
    threads=$(nproc 2>/dev/null || echo "4")
fi

# Timing
start_time=$(date +%s)

basename=$(basename "$input_vcf")
output_vcf="annotated.${basename}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting annotation process..."
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Input VCF: $input_vcf"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Output VCF: $output_vcf"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Test mode: $test_mode"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Threads: $threads"

echo ""
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finding missing FILTERs..."

# Define FILTER descriptions as Bash associative array
declare -A FILTER_DESC
FILTER_DESC["NO_HQ_GENOTYPES"]="Site has no high quality variant genotypes. No high-quality genotype (GQ>=20, DP>=10, and AB>=0.2 for heterozygotes) called for the variant. If there is one genotype at the variant site, the filter will not be applied and the variant site will pass. Allele Balance (AB) is min(AD)/DP for diploid GTs."
FILTER_DESC["ExcessHet"]="Site has excess het value larger than the threshold. Genotypes with this filter show a higher proportion of heterozygotes than expected under Hardy-Weinberg equilibrium (z-score < -4.5, phred 54.69). Suggests mapping errors or contamination."
FILTER_DESC["LowQual"]="QUALapprox is too low (lower than 60 for SNPs; lower than 69 for Indels). QUAL tells you how confident we are that there is some kind of variation at a given site."
FILTER_DESC["EXCESS_ALLELES"]="Site has an excess of alternate alleles based on the input threshold (e.g. >100 alternate alleles)."
FILTER_DESC["OUTSIDE_OF_TARGETS"]="Exome only. The site is not within the target intervals of the exome assay."

missing_header="missing_filters.txt"
> "$missing_header"

# Select (possibly temporary) VCF to operate on
vcf_for_annotation="$input_vcf"
if [[ "$test_mode" == "true" ]]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] TEST MODE: Creating subset with first 1,000 variants"
    tmp_test_file=$(mktemp --suffix=.vcf.gz)
    # Use bcftools for test mode (fast for small subsets)
    (bcftools view -h "$input_vcf"; bcftools view -H "$input_vcf" | head -n 1000) | bgzip -@"$threads" -c > "$tmp_test_file"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Indexing test subset..."
    tabix -p vcf "$tmp_test_file"
    vcf_for_annotation="$tmp_test_file"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Test subset created: $tmp_test_file"
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Extracting FILTER IDs from header..."
# Extract FILTER IDs in header (populate array)
mapfile -t header_filters < <(bcftools view -h "$vcf_for_annotation" | awk -F'[=,]' '/^##FILTER=/{print $3}')

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Found ${#header_filters[@]} FILTER headers in VCF:"
for h in "${header_filters[@]}"; do
    echo "  - $h"
done

# For each FILTER_DESC key, write a header if missing in VCF
echo ""
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Checking for missing FILTER headers..."
missing_count=0
for filter in "${!FILTER_DESC[@]}"; do
    found=0
    for h in "${header_filters[@]}"; do
        if [[ "$h" == "$filter" ]]; then
            found=1
            break
        fi
    done
    if [[ $found -eq 0 ]]; then
        echo "  [MISSING] $filter - will add to header"
        echo "##FILTER=<ID=$filter,Description=\"${FILTER_DESC[$filter]}\">" >> "$missing_header"
        ((missing_count++))
    fi
done

if [[ $missing_count -eq 0 ]]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] All required FILTER headers present"
else
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Found $missing_count missing FILTER headers"
fi

# Annotate header if missing filters present
echo ""
if [[ -s "$missing_header" ]]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Annotating VCF with missing FILTER headers..."
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 1: Creating new header..."
    
    # Step 1: Create the new header - use bcftools (fast, only reads header)
    (bcftools view -h "$vcf_for_annotation" | head -n -1;
     cat "$missing_header";
     bcftools view -h "$vcf_for_annotation" | tail -n 1) > new_header.txt
    
    header_time=$(date +%s)
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Header created (took $((header_time - start_time))s)"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 2: Appending body (all variant lines) using $threads threads..."
    
    # Get input file size for progress monitoring
    input_size=$(stat -c%s "$vcf_for_annotation" 2>/dev/null || stat -f%z "$vcf_for_annotation" 2>/dev/null)
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Input file size: $(numfmt --to=iec-i --suffix=B $input_size 2>/dev/null || echo $input_size bytes)"
    
    # Start file size monitor in background
    (
        sleep 3  # Give output file time to be created
        start_monitor=$(date +%s)
        while kill -0 $$ 2>/dev/null; do
            if [ -f "$output_vcf" ]; then
                current_size=$(stat -c%s "$output_vcf" 2>/dev/null || stat -f%z "$output_vcf" 2>/dev/null || echo 0)
                if [ "$current_size" -gt 0 ] && [ "$input_size" -gt 0 ]; then
                    pct=$((current_size * 100 / input_size))
                    elapsed=$(($(date +%s) - start_monitor))
                    if [ "$pct" -gt 0 ] && [ "$elapsed" -gt 0 ]; then
                        eta=$((elapsed * (100 - pct) / pct))
                        rate=$((current_size / elapsed))
                        printf "\r[Progress] %3d%% | %s / %s | Rate: %s/s | ETA: %dm%ds     " \
                            "$pct" \
                            "$(numfmt --to=iec-i --suffix=B $current_size 2>/dev/null || echo ${current_size}B)" \
                            "$(numfmt --to=iec-i --suffix=B $input_size 2>/dev/null || echo ${input_size}B)" \
                            "$(numfmt --to=iec-i --suffix=B $rate 2>/dev/null || echo ${rate}B)" \
                            "$((eta / 60))" "$((eta % 60))"
                    fi
                fi
            fi
            sleep 30
        done
        echo ""  # New line after progress
    ) &
    monitor_pid=$!
    
    # Step 2: Combine header + body and compress
    if command -v pigz &> /dev/null; then
        (cat new_header.txt;
         pigz -dc -p"$threads" "$vcf_for_annotation" | grep -v "^#") | bgzip -@"$threads" -c > "$output_vcf"
    else
        (cat new_header.txt;
         zcat "$vcf_for_annotation" | grep -v "^#") | bgzip -@"$threads" -c > "$output_vcf"
    fi
    
    # Stop monitor
    kill "$monitor_pid" 2>/dev/null || true
    wait "$monitor_pid" 2>/dev/null || true
    
    rm -f new_header.txt
     
    annotation_time=$(date +%s)
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Annotation complete (took $((annotation_time - start_time))s)"
else
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] No missing FILTER headers found - copying VCF..."
    cp "$vcf_for_annotation" "$output_vcf"
fi

# Index handling: only recalculate if in test mode, otherwise copy
if [[ "$test_mode" == "true" ]]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Test mode - creating new index..."
    tabix -p vcf "$output_vcf"
else
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copying original index (header-only changes)..."
    cp "${input_vcf}.tbi" "${output_vcf}.tbi"
fi

# Cleanup
if [[ "$test_mode" == "true" ]]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Cleaning up temporary test file..."
    rm -f "$tmp_test_file" "${tmp_test_file}.tbi"
fi

end_time=$(date +%s)
total_time=$((end_time - start_time))

echo ""
echo "[$(date '+%Y-%m-%d %H:%M:%S')] =========================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] COMPLETE!"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Total time: ${total_time}s"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Output: $output_vcf"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Index:  ${output_vcf}.tbi"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] =========================================="
