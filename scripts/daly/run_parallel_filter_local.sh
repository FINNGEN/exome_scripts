#!/bin/bash
# Local test script that replicates the ParallelFilter WDL task
# Usage: ./run_parallel_filter_local.sh <input.vcf.gz> [filter_expression] [cpu_count]

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <input.vcf.gz> [filter_expression] [cpu_count]"
  echo "Example: $0 input.vcf.gz 'FILTER~\"NO_HQ_GENOTYPES\"' 8"
  exit 1
fi

input_file="$1"
FILTER_EXPR="${2:-FILTER~\"NO_HQ_GENOTYPES\"}"
CHUNKS="${3:-8}"

echo "=== ParallelFilter: Starting ==="
echo "Input file: $input_file"
echo "Filter expression: $FILTER_EXPR"
echo "Target chunks: $CHUNKS"

# Check if input file exists
if [[ ! -f "$input_file" ]]; then
  echo "Error: Input file not found: $input_file"
  exit 1
fi

# Check if index exists, create if not
if [[ ! -f "$input_file.tbi" ]]; then
  echo "Creating index..."
  tabix -p vcf "$input_file"
fi

echo "Creating $CHUNKS region files..."
# Get chromosome and contig length from index
read chrom contig_len < <(bcftools index -s "$input_file" | awk '{print $1, $2}')
echo "Chromosome: $chrom, Contig length: $contig_len"

# Get first position
first_pos=$(bcftools view -H "$input_file" | head -n 1 | cut -f2)
echo "First variant position: $first_pos"

# Binary search for last variant position
echo "Starting binary search for last variant position..."
low=$first_pos; high=$contig_len
iter=0
while (( low <= high )); do
  mid=$(( (low + high) / 2 ))
  echo "  Iteration $((++iter)): Testing range $low-$high (mid=$mid)"
  if bcftools view -H -r "$chrom:$mid-$high" "$input_file" 2>/dev/null | head -n 1 | grep -q .; then
    echo "    Found variant(s), searching higher"
    low=$(( mid + 1 ))
  else
    echo "    No variants, searching lower"
    high=$(( mid - 1 ))
  fi
done

# Get exact last position from narrow window
echo "Binary search complete, narrowed to position $high"
search_start=$(( high > 10000 ? high - 10000 : first_pos ))
echo "Finding exact last position in range $search_start-$high..."
last_pos=$(bcftools view -H -r "$chrom:$search_start-$high" "$input_file" | tail -n 1 | cut -f2)
echo "Final last variant position: $last_pos"
echo "Variant range: $first_pos-$last_pos (span: $((last_pos - first_pos)) bp)"

# Split interval equally using Python linspace
echo "Splitting into $CHUNKS equal regions..."
python3 -c "import numpy as np; [open(f'region_chunk_{i:02d}','w').write(f'$chrom\t{int(s)}\t{int(e)}\n') for i,(s,e) in enumerate(zip(np.linspace($first_pos,$last_pos,$CHUNKS+1)[:-1], np.linspace($first_pos,$last_pos,$CHUNKS+1)[1:]))]"

echo "Region files created:"
for f in region_chunk_*; do
  echo "  $f: $(cat $f)"
done

# Create processing script to avoid quoting issues
cat > process_chunk.sh << 'SCRIPT_EOF'
#!/bin/bash
input_file="$1"
region_file="$2"
filter_expr="$3"
output_file="chunk${region_file}.vcf.gz"

echo "Processing $region_file"
bcftools view "$input_file" -R "$region_file" -Ou | \
  bcftools view -e "$filter_expr" -Ou | \
  bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' -Oz -o "$output_file"
echo "Completed $region_file"
SCRIPT_EOF

chmod +x process_chunk.sh

echo "Processing chunks in parallel..."
ls region_chunk_* | parallel -j $CHUNKS './process_chunk.sh '"$input_file"' {} '"'$FILTER_EXPR'"

echo "Concatenating chunks..."
bcftools concat -n -f <(ls chunk*.vcf.gz) -Oz -o filtered.vcf.gz && rm chunk*.vcf.gz region_chunk_* process_chunk.sh

echo "Indexing output..."
tabix -p vcf filtered.vcf.gz

# Count variants
variant_count=$(bcftools view -H filtered.vcf.gz | wc -l)
echo "Total variants in filtered file: $variant_count"

echo "=== ParallelFilter: Done! Output file: filtered.vcf.gz ==="
