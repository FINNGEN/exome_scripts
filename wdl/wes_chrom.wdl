version 1.0

workflow wes_chrom {
  input {
    File vcf_list
    String genotype_filter
    String variant_filter
    Int cpu_count
    Int? test_sample_count
    String output_prefix
   }

  Array[File] vcf_files = read_lines(vcf_list)
  
  scatter (vcf in vcf_files) {
    if (defined(test_sample_count)) {
      call SubsetSamples {
        input:
          input_vcf = vcf,
          sample_count = select_first([test_sample_count])
      }
    }

    File vcf_to_filter = select_first([SubsetSamples.subset_vcf, vcf])

    call ComputeStats as OriginalStats {
      input:
        input_vcf = vcf_to_filter,
        cpu_count = cpu_count
    }
    
    call ParallelFilterByRegion {
      input:
        input_vcf = vcf_to_filter,
        positions = OriginalStats.positions,
        genotype_filter = genotype_filter,
        variant_filter = variant_filter,
        cpu_count = cpu_count
    }
    
    call ComputeStats as FilteredStats {
      input:
        input_vcf = ParallelFilterByRegion.filtered_vcf,
        cpu_count = cpu_count
    }
    
    call ValidateFiltering {
      input:
        original_sample_vcf = OriginalStats.sample_vcf,
        original_sample_vcf_tbi = OriginalStats.sample_vcf_tbi,
        original_stats = OriginalStats.stats,
        filtered_sample_vcf = FilteredStats.sample_vcf,
        filtered_sample_vcf_tbi = FilteredStats.sample_vcf_tbi,
        filtered_stats = FilteredStats.stats,
        genotype_filter = genotype_filter,
        variant_filter = variant_filter
    }
  }

  call SortAndMerge {
    input:
    vcf_files = ParallelFilterByRegion.filtered_vcf,
    vcf_tbi_files = ParallelFilterByRegion.filtered_vcf_tbi,
    cpu_count = cpu_count,
    output_prefix = output_prefix
  }

  call SummaryStats {
    input:
      vcf_file_names = vcf_files,
      original_stats = OriginalStats.stats,
      filtered_stats = FilteredStats.stats
  }

  output {
    Array[File] filtered_vcfs = ParallelFilterByRegion.filtered_vcf
    Array[File] filtered_vcf_tbis = ParallelFilterByRegion.filtered_vcf_tbi
    Array[File] original_stats = OriginalStats.stats
    Array[File] filtered_stats = FilteredStats.stats
    Array[File] validation_reports = ValidateFiltering.report
    File merged_vcf = SortAndMerge.merged_vcf
    File merged_vcf_tbi = SortAndMerge.merged_vcf_tbi
    File summary_table = SummaryStats.summary
  }
}

task ComputeStats {
  input {
    File input_vcf
    Int cpu_count = 8
  }

  File input_vcf_tbi = input_vcf + ".tbi"
  Int disk_size = ceil(size(input_vcf, 'GB')) + 10

  command <<<
  set -euo

  echo "=== Computing statistics and creating sample ==="
  
  CHUNKS=~{cpu_count}
  
  # Touch index to ensure it's localized
  touch ~{input_vcf_tbi}

  # Get chromosome from first variant
  chrom=$(bcftools query -f '%CHROM\n' "~{input_vcf}" | head -n 1)
  echo "Chromosome: $chrom"

  # Get chromosome and contig length from index
  read chrom_idx contig_len < <(bcftools index -s "~{input_vcf}" | awk '{print $1, $2}')
  echo "Contig length: $contig_len"

  # Get first position
  first_pos=$(bcftools view -H "~{input_vcf}" | head -n 1 | cut -f2)
  echo "First variant position: $first_pos"

  # Binary search for last variant position
  echo "Starting binary search for last variant..."
  low=$first_pos
  high=$contig_len
  while (( low <= high )); do
    mid=$(( (low + high) / 2 ))
    if bcftools view -H -r "$chrom:$mid-$high" "~{input_vcf}" 2>/dev/null | head -n 1 | grep -q .; then
      low=$(( mid + 1 ))
    else
      high=$(( mid - 1 ))
    fi
  done

  # Get exact last position from narrow window
  search_start=$(( high > 10000 ? high - 10000 : first_pos ))
  last_pos=$(bcftools view -H -r "$chrom:$search_start-$high" "~{input_vcf}" | tail -n 1 | cut -f2)
  echo "Last variant position: $last_pos"
  echo ""

  # Split interval equally using Python linspace
  echo "Splitting into $CHUNKS equal regions for parallel processing..."
  python3 -c "import numpy as np; [open(f'region_chunk_{i:02d}','w').write(f'$chrom\t{int(s)}\t{int(e)}\n') for i,(s,e) in enumerate(zip(np.linspace($first_pos,$last_pos,$CHUNKS+1)[:-1], np.linspace($first_pos,$last_pos,$CHUNKS+1)[1:]))]"

  # Create processing script
  cat > extract_chunk.sh << 'SCRIPT_EOF'
  #!/bin/bash
  input_file="$1"
  region_file="$2"
  
  bcftools query -f '%POS\n' "$input_file" -R "$region_file" > "${region_file}.positions"
  SCRIPT_EOF
  chmod +x extract_chunk.sh

  echo "Extracting positions in parallel..."
  ls region_chunk_* | sort -V | parallel -j $CHUNKS './extract_chunk.sh "~{input_vcf}" {}'

  echo "Concatenating and sorting position files..."
  cat region_chunk_*.positions | sort -n -u > positions.txt
  
  # Count total variants
  variant_count=$(wc -l < positions.txt)
  echo "Total variants: $variant_count"
  
  # Cleanup intermediate files
  rm -f region_chunk_*.positions extract_chunk.sh region_chunk_*
  
  # Create a small sample VCF (100 variants) for validation
  echo "Creating sample VCF for validation (100 variants)..."
  bcftools view -h "~{input_vcf}" | bgzip -c > sample.vcf.gz
  bcftools view -H "~{input_vcf}" | head -n 100 | bgzip -c >> sample.vcf.gz
  tabix -p vcf sample.vcf.gz
  echo "Sample VCF created: $(bcftools view -H sample.vcf.gz | wc -l) variants"
  
  # Create stats file
  echo -e "variant_count\t$variant_count" > stats.txt
  
  echo "=== Complete ==="
  >>>

  output {
    File positions = "positions.txt"
    File sample_vcf = "sample.vcf.gz"
    File sample_vcf_tbi = "sample.vcf.gz.tbi"
    File stats = "stats.txt"
  }

  runtime {
    memory: "4G"
    disks: "local-disk ~{disk_size} HDD"
    cpu: cpu_count
    preemptible: 1
  }
}

task ParallelFilterByRegion {
  input {
    File input_vcf
    File positions
    String genotype_filter
    String variant_filter
    Int cpu_count
  }

  File input_vcf_tbi = input_vcf + ".tbi"
  String base_name = basename(basename(basename(input_vcf, ".vcf.gz"), ".vcf.bgz"), ".bcf")
  Int disk_size = ceil(size(input_vcf,'GB')*3) + 20
  Int vcf_size_gb = ceil(size(input_vcf, 'GB'))
  Int memory_gb = if vcf_size_gb <= 8 then 8
                  else if vcf_size_gb <= 16 then 16
                  else if vcf_size_gb <= 32 then 32
                  else if vcf_size_gb <= 64 then 64
                  else 128
  
  command <<<
  set -euo
  
  input_file="~{input_vcf}"
  CHUNKS=~{cpu_count}
  touch ~{input_vcf_tbi}

  echo "=== Parallel Filter by Region ==="
  echo "Input: $input_file"
  echo "Genotype filter: ~{genotype_filter}"
  echo "Variant filter: ~{variant_filter}"
  echo "CPU cores: $CHUNKS"
  echo ""
  
  # Get chromosome from first variant
  chrom=$(bcftools query -f '%CHROM\n' "$input_file" | head -n 1)
  echo "Chromosome: $chrom"
  echo ""
  
  # Split positions into equal chunks locally
  echo "Splitting positions into $CHUNKS equal chunks..."
  split -n l/$CHUNKS -d ~{positions} pos_chunk_tmp_
  
  # Create properly formatted position and region files
  chunk_num=0
  for tmp_file in pos_chunk_tmp_*; do
    variant_count=$(wc -l < "$tmp_file")
    first_pos=$(head -n 1 "$tmp_file")
    last_pos=$(tail -n 1 "$tmp_file")
    
    chunk_id=$(printf '%02d' $chunk_num)
    pos_name="pos_chunk_${chunk_id}"
    region_name="region_chunk_${chunk_id}"
    
    # Create position file with CHROM POS format (required by -T)
    awk -v chr="$chrom" '{print chr "\t" $1}' "$tmp_file" > "$pos_name"
    
    # Create region file for efficient seeking
    printf "%s:%d-%d\n" "$chrom" "$first_pos" "$last_pos" > "$region_name"
    
    echo "  chunk_${chunk_id}: $chrom:$first_pos-$last_pos ($variant_count variants)"
    chunk_num=$((chunk_num + 1))
  done
  rm -f pos_chunk_tmp_*
  echo ""

  # Create commands file with each command fully self-contained
  echo "Creating commands for parallel processing..."
  rm -f commands.txt
  chunk_num=0
  for pos_file in pos_chunk_*; do
    chunk_id=$(printf '%02d' $chunk_num)
    region=$(cat "region_chunk_${chunk_id}")
    cat >> commands.txt << EOF
echo "Processing: chunk_${chunk_id} (region: $region)" && bcftools view "$input_file" -r "$region" -T "pos_chunk_${chunk_id}" -Ou | bcftools +setGT -Ou -- -t q -n . -i '~{genotype_filter}' | bcftools +fill-tags -Ou -- -t AC | bcftools view -i '~{variant_filter}' -Ou | bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' -Oz -o "chunk_${chunk_id}.vcf.gz" && echo "  ✓ Done: chunk_${chunk_id}"
EOF
    chunk_num=$((chunk_num + 1))
  done
  
  echo "Generated commands (first 3):"
  head -3 commands.txt
  echo ""

  echo "Processing chunks in parallel..."
  cat commands.txt | parallel -j $CHUNKS

  echo ""
  echo "Concatenating chunks..."
  ls chunk_*.vcf.gz | sort -V > chunk_list.txt
  bcftools concat -n -f chunk_list.txt -Oz -o ~{base_name}.filtered.vcf.gz
  
  echo "Indexing output..."
  tabix -p vcf ~{base_name}.filtered.vcf.gz
  
  # Cleanup
  echo "Cleaning up temporary files..."
  rm -f chunk_*.vcf.gz chunk_list.txt pos_chunk_* region_chunk_* commands.txt
  
  echo "=== Complete! ==="
  >>>

  output {
    File filtered_vcf = "~{base_name}.filtered.vcf.gz"
    File filtered_vcf_tbi = "~{base_name}.filtered.vcf.gz.tbi"
  }

  runtime {
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_size} HDD"
    cpu: cpu_count
    preemptible: 1
  }
}

task SortAndMerge {
  input {
    Array[File] vcf_files
    Array[File] vcf_tbi_files
    Int cpu_count
    String output_prefix
  }

  Int disk_size = ceil(size(vcf_files, 'GB') * 2) + 50
  String out_vcf = output_prefix + ".merged.vcf.gz"
  String out_tbi = output_prefix + ".merged.vcf.gz.tbi"
  
  command <<<
  set -euo

  echo "=== Sorting and Merging VCF Files ==="
  echo "Number of input files: ~{length(vcf_files)}"
  
  # Create file list and sort by version
  sort -V ~{write_lines(vcf_files)} > vcf_list.txt

  # Touch all indices
  while read -r tbi; do
      touch "$tbi"
  done < <(cat ~{write_lines(vcf_tbi_files)})

  echo "Sorting VCF files by chromosome and position..."
  bcftools concat -f vcf_list.txt -Oz -o ~{out_vcf} --threads ~{cpu_count}

  echo "Indexing merged VCF..."
  tabix -p vcf ~{out_vcf}
  
  echo "=== Merge Complete! ==="
  >>>
  output {
    File merged_vcf = out_vcf
    File merged_vcf_tbi = out_tbi
  }

  runtime {
    memory: "16G"
    disks: "local-disk ~{disk_size} HDD"
    cpu: cpu_count
    preemptible: 1
  }
}

task ValidateFiltering {
  input {
    File original_sample_vcf
    File original_sample_vcf_tbi
    File original_stats
    File filtered_sample_vcf
    File filtered_sample_vcf_tbi
    File filtered_stats
    String genotype_filter
    String variant_filter
  }

  Int disk_size = 10

  command <<<
  set -euo

  echo "=== Validating VCF Filtering ===" > report.txt
  echo "Sample-based validation (using up to 10k variants)" >> report.txt
  echo "" >> report.txt

  # Touch indices
  touch ~{original_sample_vcf_tbi}
  touch ~{filtered_sample_vcf_tbi}

  # Read variant counts
  orig_count=$(grep "variant_count" ~{original_stats} | cut -f2)
  filt_count=$(grep "variant_count" ~{filtered_stats} | cut -f2)
  
  echo "Original variant count: $orig_count" >> report.txt
  echo "Filtered variant count: $filt_count" >> report.txt
  
  if [[ "$orig_count" -gt 0 ]]; then
    pct_kept=$(awk "BEGIN {printf \"%.1f\", ($filt_count / $orig_count) * 100}")
    pct_dropped=$(awk "BEGIN {printf \"%.1f\", (($orig_count - $filt_count) / $orig_count) * 100}")
  else
    pct_kept="0.0"
    pct_dropped="0.0"
  fi
  
  echo "Kept: ${pct_kept}%" >> report.txt
  echo "Dropped: ${pct_dropped}%" >> report.txt
  echo "" >> report.txt

  # Test 1: Check genotypes were set to missing
  echo "Test 1: Checking genotypes were set to missing (~{genotype_filter})" >> report.txt
  echo "--------------------------------------------------------------" >> report.txt

  found=0
  bcftools query -i 'GT!="mis"' -f '%CHROM\t%POS[\t%SAMPLE\t%GT\t%DP\t%GQ]\n' "~{original_sample_vcf}" 2>/dev/null | \
    awk '($4 < 10 || $5 < 20)' | head -100 | \
    while read chrom pos sample gt dp gq; do
      if bcftools view -H -r "$chrom:$pos" "~{filtered_sample_vcf}" 2>/dev/null | grep -q .; then
        after_gt=$(bcftools query -s "$sample" -r "$chrom:$pos" -f '[\t%GT]\n' "~{filtered_sample_vcf}" 2>/dev/null | tr -d '\t' | tr -d ' ')
        
        if [[ "$after_gt" == "./." ]] || [[ "$after_gt" == "." ]]; then
          echo "✓ PASS: $chrom:$pos sample=$sample GT: $gt → $after_gt" >> report.txt
        else
          echo "✗ FAIL: $chrom:$pos sample=$sample GT: $gt → $after_gt (not missing)" >> report.txt
        fi
        found=1
        break
      fi
    done

  if [[ $found -eq 0 ]]; then
    echo "⚠ NOTE: All low-quality genotypes resulted in variants being removed (AC=0)" >> report.txt
  fi

  echo "" >> report.txt

  # Test 2: Check AC field exists
  echo "Test 2: Checking AC field is present in filtered VCF" >> report.txt
  echo "-----------------------------------------------------" >> report.txt

  first_var=$(bcftools view -H "~{filtered_sample_vcf}" | head -1 | awk '{print $1":"$2}')

  if [[ -n "$first_var" ]]; then
    ac_value=$(bcftools query -r "$first_var" -f '%AC\n' "~{filtered_sample_vcf}" 2>/dev/null)
    
    if [[ -n "$ac_value" ]]; then
      echo "✓ PASS: AC field present (example: $first_var AC=$ac_value)" >> report.txt
    else
      echo "✗ FAIL: AC field missing in filtered VCF" >> report.txt
    fi
  else
    echo "⚠ NOTE: No variants in filtered file" >> report.txt
  fi

  echo "" >> report.txt

  # Test 3: Check variant IDs
  echo "Test 3: Checking variant IDs (CHROM_POS_REF_ALT format)" >> report.txt
  echo "--------------------------------------------------------" >> report.txt

  sample_id=$(bcftools view -H "~{filtered_sample_vcf}" 2>/dev/null | head -1 | awk '{print $3}')
  if [[ "$sample_id" =~ ^[^_]+_[0-9]+_.+_.+$ ]]; then
    echo "✓ PASS: IDs formatted as CHROM_POS_REF_ALT (example: $sample_id)" >> report.txt
  else
    echo "✗ FAIL: ID not in expected format (got: $sample_id)" >> report.txt
  fi

  echo "" >> report.txt
  echo "=== Validation Complete ===" >> report.txt
  >>>

  output {
    File report = "report.txt"
  }

  runtime {
    memory: "4G"
    disks: "local-disk ~{disk_size} HDD"
    cpu: 1
    preemptible: 1
  }
}

task SubsetSamples {
  input {
    File input_vcf
    Int sample_count
  }

  String base_name = basename(basename(basename(input_vcf, ".vcf.gz"), ".vcf.bgz"), ".bcf")
  String output_vcf = base_name + ".subset_" + sample_count + "samples.vcf.gz"
  String output_tbi = base_name + ".subset_" + sample_count + "samples.vcf.gz.tbi"
  Int disk_size = ceil(size(input_vcf, 'GB') * 2) + 20

  command <<<
  set -euo
    
  THREADS=$(nproc)
  NCOLS=$((9 + ~{sample_count}))
    
  echo "=== Subsetting VCF to first ~{sample_count} samples (columns 1-$NCOLS) ==="
  
  # Extract metadata lines (##) - keep intact
  echo "Processing metadata..."
  bcftools view -h "~{input_vcf}" | grep "^##" | bgzip -@ $THREADS -c > ~{output_vcf}
  
  # Extract and cut column header line (#CHROM)
  echo "Processing column header..."
  bcftools view -h "~{input_vcf}" | grep "^#CHROM" | cut -f 1-$NCOLS | bgzip -@ $THREADS -c >> ~{output_vcf}
  
  # Extract and cut body (skip all header lines starting with #)
  echo "Processing body..."
  zcat "~{input_vcf}" | grep -v "^#" | cut -f 1-$NCOLS | bgzip -@ $THREADS -c >> ~{output_vcf}
  
  echo "Indexing..."
  tabix -p vcf ~{output_vcf}
  echo "=== Complete ==="
  >>>

  output {
    File subset_vcf = output_vcf
    File subset_vcf_tbi = output_tbi
  }

  runtime {
    memory: "8G"
    disks: "local-disk ~{disk_size} HDD"
    cpu: 4
    preemptible: 1
  }
}

task SummaryStats {
  input {
    Array[String] vcf_file_names
    Array[File] original_stats
    Array[File] filtered_stats
  }

  Int disk_size = 10

  command <<<
  set -euo

  echo "=== Creating Summary Statistics ==="
  
  # Create header
  echo -e "chromosome\toriginal_variants\tfiltered_variants\tpercent_dropped" > summary.tsv
  
  # Process each VCF file
  idx=0
  while IFS= read -r vcf_file; do
    # Extract chromosome from filename (e.g., chr1.vcf.gz -> chr1)
    filename=$(basename "$vcf_file")
    chrom=$(echo "$filename" | sed -E 's/.*_(chr[0-9XY]+).*\.vcf\.(gz|bgz)$/\1/' | sed -E 's/^.*chr([0-9XY]+).*$/chr\1/' | sed -E 's/.*[^a-zA-Z](chr[0-9XY]+)[^a-zA-Z].*$/\1/' | head -1)
    
    # Fallback: try to extract just the chromosome number/letter
    if [[ -z "$chrom" ]] || [[ "$chrom" == "$filename" ]]; then
      chrom=$(echo "$filename" | grep -oE '(chr)?[0-9XY]+' | head -1)
    fi
    
    # If still empty, use filename
    if [[ -z "$chrom" ]]; then
      chrom=$(echo "$filename" | sed 's/\..*//')
    fi
    
    # Read original and filtered counts from stats files (not VCFs)
    orig_stats_file=$(echo "~{sep=' ' original_stats}" | cut -d' ' -f$((idx+1)))
    filt_stats_file=$(echo "~{sep=' ' filtered_stats}" | cut -d' ' -f$((idx+1)))
    
    orig_count=$(grep "variant_count" "$orig_stats_file" | cut -f2)
    filt_count=$(grep "variant_count" "$filt_stats_file" | cut -f2)
    
    # Calculate percent dropped
    if [[ "$orig_count" -gt 0 ]]; then
      pct_dropped=$(awk "BEGIN {printf \"%.2f\", (($orig_count - $filt_count) / $orig_count) * 100}")
    else
      pct_dropped="0.00"
    fi
    
    echo -e "${chrom}\t${orig_count}\t${filt_count}\t${pct_dropped}" >> summary.tsv
    idx=$((idx + 1))
  done < <(cat << 'EOF'
~{sep='\n' vcf_file_names}
EOF
)
  
  # Add totals row
  echo "" >> summary.tsv
  total_orig=$(grep "variant_count" ~{sep=' ' original_stats} | cut -f2 | awk '{sum+=$1} END {print sum}')
  total_filt=$(grep "variant_count" ~{sep=' ' filtered_stats} | cut -f2 | awk '{sum+=$1} END {print sum}')
  total_pct=$(awk "BEGIN {printf \"%.2f\", (($total_orig - $total_filt) / $total_orig) * 100}")
  
  echo -e "TOTAL\t${total_orig}\t${total_filt}\t${total_pct}" >> summary.tsv
  
  echo ""
  echo "Summary Table:"
  column -t summary.tsv
  
  echo "=== Complete ==="
  >>>

  output {
    File summary = "summary.tsv"
  }

  runtime {
    memory: "2G"
    disks: "local-disk ~{disk_size} HDD"
    cpu: 1
    preemptible: 1
  }
}

