version 1.0

workflow daly_qc {
  input {
    File vcf_list
    String filter_expression = 'FILTER~"NO_HQ_GENOTYPES"'
    Int cpu_count = 8
  }
  Array[File] vcf_files = read_lines(vcf_list)
  scatter (vcf in vcf_files) {
    call ComputeStats as OriginalStats {
      input:
        input_vcf = vcf
    }
    
    call AnnotateHeaders {
      input:
        input_vcf = vcf,
        cpu_count = cpu_count
    }
    
    call ParallelFilter {
      input:
        input_vcf = AnnotateHeaders.annotated_vcf,
        input_vcf_tbi = AnnotateHeaders.annotated_vcf_tbi,
        filter_expression = filter_expression,
        cpu_count = cpu_count
    }
    
    call ComputeStats as FilteredStats {
      input:
        input_vcf = ParallelFilter.filtered_vcf
    }
    
    call ValidateFiltering {
      input:
        original_sample_vcf = OriginalStats.sample_vcf,
        original_sample_vcf_tbi = OriginalStats.sample_vcf_tbi,
        original_stats = OriginalStats.stats,
        filtered_sample_vcf = FilteredStats.sample_vcf,
        filtered_sample_vcf_tbi = FilteredStats.sample_vcf_tbi,
        filtered_stats = FilteredStats.stats,
        filter_expression = filter_expression
    }
  }

  call SummaryStats {
    input:
      vcf_file_names = vcf_files,
      original_stats = OriginalStats.stats,
      filtered_stats = FilteredStats.stats
  }

  call SortAndMerge {
    input:
      vcf_files = ParallelFilter.filtered_vcf,
      vcf_tbi_files = ParallelFilter.filtered_vcf_tbi,
      summary_report = SummaryStats.report,
      root_name = SummaryStats.root_name,
      cpu_count = cpu_count
  }

  output {
    Array[File] filtered_vcfs = ParallelFilter.filtered_vcf
    Array[File] filtered_vcf_tbis = ParallelFilter.filtered_vcf_tbi
    Array[File] original_stats = OriginalStats.stats
    Array[File] filtered_stats = FilteredStats.stats
    Array[File] validation_reports = ValidateFiltering.report
    File merged_vcf = SortAndMerge.merged_vcf
    File merged_vcf_tbi = SortAndMerge.merged_vcf_tbi
    File report = SortAndMerge.report
  }
}

task AnnotateHeaders {
  input {
    File input_vcf
    Int cpu_count = 4
  }

  Int disk_size = ceil(size(input_vcf,'GB')*2) + 20
  command <<<
  
  basename=$(basename "~{input_vcf}")
  output_vcf="annotated.${basename}"

  echo "Finding missing FILTERs..."

  # Define FILTER descriptions as Bash associative array
  declare -A FILTER_DESC
  FILTER_DESC["NO_HQ_GENOTYPES"]="Site has no high quality variant genotypes. No high-quality genotype (GQ>=20, DP>=10, and AB>=0.2 for heterozygotes) called for the variant. If there is one genotype at the variant site, the filter will not be applied and the variant site will pass. Allele Balance (AB) is min(AD)/DP for diploid GTs."
  FILTER_DESC["ExcessHet"]="Site has excess het value larger than the threshold. Genotypes with this filter show a higher proportion of heterozygotes than expected under Hardy-Weinberg equilibrium (z-score < -4.5, phred 54.69). Suggests mapping errors or contamination."
  FILTER_DESC["LowQual"]="QUALapprox is too low (lower than 60 for SNPs; lower than 69 for Indels). QUAL tells you how confident we are that there is some kind of variation at a given site."
  FILTER_DESC["EXCESS_ALLELES"]="Site has an excess of alternate alleles based on the input threshold (e.g. >100 alternate alleles)."
  FILTER_DESC["OUTSIDE_OF_TARGETS"]="Exome only. The site is not within the target intervals of the exome assay."

  missing_header="missing_filters.txt"
  > "$missing_header"

  # Extract FILTER IDs in header (populate array)
  mapfile -t header_filters < <(bcftools view -h "~{input_vcf}" | awk -F'[=,]' '/^##FILTER=/{print $3}')
  
  # For each FILTER_DESC key, write a header if missing in VCF
  for filter in "${!FILTER_DESC[@]}"; do
      found=0
      for h in "${header_filters[@]}"; do
          if [[ "$h" == "$filter" ]]; then
              found=1
              break
          fi
      done
      if [[ $found -eq 0 ]]; then
          echo "##FILTER=<ID=$filter,Description=\"${FILTER_DESC[$filter]}\">" >> "$missing_header"
      fi
  done

  # Annotate header if missing filters present
  if [[ -s "$missing_header" ]]; then
      echo "Annotating with missing FILTER headers..."
      # Step 1: Create new header with bcftools (fast, only reads header)
      (bcftools view -h "~{input_vcf}" | head -n -1;
       cat "$missing_header";
       bcftools view -h "~{input_vcf}" | tail -n 1) > new_header.txt
      
      echo "Appending body with parallel compression..."
      
      # Get input file size for progress monitoring
      input_size=$(stat -c%s "~{input_vcf}" 2>/dev/null || stat -f%z "~{input_vcf}" 2>/dev/null)
      echo "Input file size: $(numfmt --to=iec-i --suffix=B $input_size 2>/dev/null || echo $input_size bytes)"
      
      # Start file size monitor in background
      (
          sleep 3
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
          echo ""
      ) &
      monitor_pid=$!
      
      # Step 2: Append body with zcat (fast, no VCF parsing)
      if command -v pigz &> /dev/null; then
          echo "Using pigz for parallel decompression"
          (cat new_header.txt;
           pigz -dc -p~{cpu_count} "~{input_vcf}" | grep -v "^#") | bgzip -@~{cpu_count} -c > "$output_vcf"
      else
          (cat new_header.txt;
           zcat "~{input_vcf}" | grep -v "^#") | bgzip -@~{cpu_count} -c > "$output_vcf"
      fi
      
      # Stop monitor
      kill "$monitor_pid" 2>/dev/null || true
      wait "$monitor_pid" 2>/dev/null || true
      
      rm -f new_header.txt
      
      # Index the modified VCF
      echo "Indexing modified VCF..."
      tabix -p vcf "$output_vcf"
  else
      echo "No missing headers - copying VCF and index..."
      cp "~{input_vcf}" "$output_vcf"
      # Copy the index instead of recalculating (header-only changes don't affect index)
      cp "~{input_vcf}.tbi" "${output_vcf}.tbi"
  fi
  
  >>>

  output {
    File annotated_vcf = "annotated.${basename(input_vcf)}"
    File annotated_vcf_tbi = "annotated.${basename(input_vcf)}.tbi"
  }

  runtime {
    memory: "4G"
    disks: "local-disk ~{disk_size} HDD"
    cpu: cpu_count
    preemptible: 1
  }
}

task ParallelFilter {
  input {
    File input_vcf
    File input_vcf_tbi
    String filter_expression
    Int cpu_count
  }

  Int disk_size = ceil(size(input_vcf,'GB')*3) + 20
  
  command <<<
  
  input_file="~{input_vcf}"
  CHUNKS=~{cpu_count}
  FILTER_EXPR='~{filter_expression}'
  touch ~{input_vcf_tbi} # Ensure index is present for bcftools indexing

  echo "Creating $CHUNKS region files..."
  # Get chromosome and contig length from index
  read chrom contig_len < <(bcftools index -s "$input_file" | awk '{print $1, $2}')
  # Get first position
  first_pos=$(bcftools view -H "$input_file" | head -n 1 | cut -f2)
  # Binary search for last variant position
  low=$first_pos; high=$contig_len
  while (( low <= high )); do
    mid=$(( (low + high) / 2 ))
    if bcftools view -H -r "$chrom:$mid-$high" "$input_file" 2>/dev/null | head -n 1 | grep -q .; then
      low=$(( mid + 1 ))
    else
      high=$(( mid - 1 ))
    fi
  done
  # Get exact last position from narrow window
  search_start=$(( high > 10000 ? high - 10000 : first_pos ))
  last_pos=$(bcftools view -H -r "$chrom:$search_start-$high" "$input_file" | tail -n 1 | cut -f2)
  echo "Final last variant position: $last_pos"

  # Split interval equally using Python linspace
  python3 -c "import numpy as np; [open(f'region_chunk_{i:02d}','w').write(f'$chrom\t{int(s)}\t{int(e)}\n') for i,(s,e) in enumerate(zip(np.linspace($first_pos,$last_pos,$CHUNKS+1)[:-1], np.linspace($first_pos,$last_pos,$CHUNKS+1)[1:]))]"

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
  bcftools concat -n -f <(ls chunk*.vcf.gz) -Oz -o filtered.vcf.gz && rm chunk*.vcf.gz region_chunk_*

  echo "Indexing output..."
  tabix -p vcf filtered.vcf.gz
  
  # Count variants
  variant_count=$(bcftools view -H filtered.vcf.gz | wc -l)
  echo "Total variants in filtered file: $variant_count"
  echo "$variant_count" > variant_count.txt
  
  echo "Done! Output file: filtered.vcf.gz"
  >>>

  output {
    File filtered_vcf = "filtered.vcf.gz"
    File filtered_vcf_tbi = "filtered.vcf.gz.tbi"
    Int variant_count = read_int("variant_count.txt")
  }

  runtime {
    memory: "8G"
    disks: "local-disk ~{disk_size} HDD"
    cpu: cpu_count
  }
}

task SortAndMerge {
  input {
    Array[File] vcf_files
    Array[File] vcf_tbi_files
    File summary_report
    String root_name
    Int cpu_count = 8
  }

  Int disk_size = ceil(size(vcf_files,'GB')*2) + 50

  command <<<

  # Create file list
  cat ~{write_lines(vcf_files)} > unsorted_vcf_list.txt

  # Sort by filename (handles mixed cached/fresh GCS paths)
  awk -F'/' '{print $NF"\t"$0}' unsorted_vcf_list.txt | sort -V | cut -f2- > sorted_vcf_list.txt

  echo "Sorted VCF files:"
  cat sorted_vcf_list.txt
  echo "Output filename: ~{root_name}.QC_ANNOTATED.vcf.gz"

  echo "Concatenating sorted chromosomes..."
  bcftools concat -f sorted_vcf_list.txt -Oz -o ~{root_name}.QC_ANNOTATED.vcf.gz

  echo "Indexing merged VCF..."
  tabix -p vcf ~{root_name}.QC_ANNOTATED.vcf.gz

  cp ~{summary_report} ~{root_name}.QC_ANNOTATED.report.txt

  echo "Done! Output file: ~{root_name}.QC_ANNOTATED.vcf.gz"
  >>>

  output {
    File merged_vcf = "~{root_name}.QC_ANNOTATED.vcf.gz"
    File merged_vcf_tbi = "~{root_name}.QC_ANNOTATED.vcf.gz.tbi"
    File report = "~{root_name}.QC_ANNOTATED.report.txt"
  }

  runtime {
    memory: "8G"
    disks: "local-disk ~{disk_size} HDD"
    cpu: cpu_count
    preemptible: 1
  }
}

task ComputeStats {
  input {
    File input_vcf
  }

  File input_vcf_tbi = input_vcf + ".tbi"
  Int disk_size = ceil(size(input_vcf, 'GB')) + 10

  command <<<
  set -euo

  echo "=== Computing variant count ==="
  
  # Touch index to ensure it's localized
  touch ~{input_vcf_tbi}
  
  # Get chromosome from first variant
  chrom=$(bcftools query -f '%CHROM\n' "~{input_vcf}" | head -n 1)
  echo "Chromosome: $chrom"

  # Count total variants
  variant_count=$(bcftools view -H "~{input_vcf}" | wc -l)
  echo "Total variants: $variant_count"
  
  # Create a small sample VCF (100 variants) for chromosome validation
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
    File stats = "stats.txt"
    File sample_vcf = "sample.vcf.gz"
    File sample_vcf_tbi = "sample.vcf.gz.tbi"
  }

  runtime {
    memory: "2G"
    disks: "local-disk ~{disk_size} HDD"
    cpu: 1
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

  # Derive root name from first VCF: strip _chrN.vcf.gz (or .vcf.bgz/.bcf)
  first_vcf=$(head -1 ~{write_lines(vcf_file_names)})
  filename=$(basename "$first_vcf")
  root=$(echo "$filename" | sed -E 's/_chr[0-9XYxy]+\.(vcf\.gz|vcf\.bgz|bcf)$//')
  echo "$root" > root_name.txt
  echo "Root name: $root"

  # Create header
  echo -e "chromosome\toriginal_variants\tfiltered_variants\tpercent_dropped" > summary.report.txt
  
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
    
    echo -e "${chrom}\t${orig_count}\t${filt_count}\t${pct_dropped}" >> summary.report.txt
    idx=$((idx + 1))
  done < <(cat << 'EOF'
~{sep='\n' vcf_file_names}
EOF
)
  
  # Add totals row
  echo "" >> summary.report.txt
  total_orig=$(grep "variant_count" ~{sep=' ' original_stats} | cut -f2 | awk '{sum+=$1} END {print sum}')
  total_filt=$(grep "variant_count" ~{sep=' ' filtered_stats} | cut -f2 | awk '{sum+=$1} END {print sum}')
  total_pct=$(awk "BEGIN {printf \"%.2f\", (($total_orig - $total_filt) / $total_orig) * 100}")
  
  echo -e "TOTAL\t${total_orig}\t${total_filt}\t${total_pct}" >> summary.report.txt
  
  echo ""
  echo "Summary Table:"
  column -t summary.report.txt
  
  echo "=== Complete ==="
  >>>

  output {
    String root_name = read_string("root_name.txt")
    File report = "summary.report.txt"
  }

  runtime {
    memory: "2G"
    disks: "local-disk ~{disk_size} HDD"
    cpu: 1
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
    String filter_expression
  }

  Int disk_size = 10

  command <<<
  set -euo

  echo "=== Validating VCF Filtering ===" > report.txt
  echo "Sample-based validation (using up to 100 variants)" >> report.txt
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

  # Test 1: Check variants were filtered out based on FILTER expression
  echo "Test 1: Checking variants matching filter expression were removed" >> report.txt
  echo "Filter expression: ~{filter_expression}" >> report.txt
  echo "-------------------------------------------------------------------" >> report.txt

  # Check if any variants in filtered file match the filter expression
  filtered_matching=$(bcftools view -H -e '~{filter_expression}' "~{filtered_sample_vcf}" 2>/dev/null | wc -l || echo "0")
  filtered_total=$(bcftools view -H "~{filtered_sample_vcf}" 2>/dev/null | wc -l || echo "0")
  
  if [[ "$filtered_matching" -eq "$filtered_total" ]]; then
    echo "✓ PASS: No variants matching filter expression found in filtered file" >> report.txt
  elif [[ "$filtered_matching" -gt 0 ]]; then
    echo "✗ FAIL: Found $filtered_matching variants matching filter expression (should be 0)" >> report.txt
  else
    echo "⚠ NOTE: Unable to check filter expression" >> report.txt
  fi

  echo "" >> report.txt

  # Test 2: Check variant IDs
  echo "Test 2: Checking variant IDs (CHROM_POS_REF_ALT format)" >> report.txt
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
    memory: "2G"
    disks: "local-disk ~{disk_size} HDD"
    cpu: 1
    preemptible: 1
  }
}
