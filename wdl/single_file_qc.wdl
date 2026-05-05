version 1.0

workflow single_file_qc {
  input {
    Array[File] vcf_files
    String genotype_filter
    String variant_filter
    Int cpu_count = 8
  }

  scatter (vcf in vcf_files) {
    call FilterByChromosome {
      input:
        input_vcf = vcf,
        genotype_filter = genotype_filter,
        variant_filter = variant_filter,
        cpu_count = cpu_count
    }

    call ValidateFiltering {
      input:
        original_vcf = vcf,
        filtered_vcf = FilterByChromosome.filtered_vcf,
        filtered_vcf_tbi = FilterByChromosome.filtered_vcf_tbi,
        genotype_filter = genotype_filter
    }
  }

  output {
    Array[File] filtered_vcfs = FilterByChromosome.filtered_vcf
    Array[File] filtered_vcf_tbis = FilterByChromosome.filtered_vcf_tbi
    Array[File] validation_reports = ValidateFiltering.report
  }
}

task FilterByChromosome {
  input {
    File input_vcf
    String genotype_filter
    String variant_filter
    Int cpu_count
  }

  File input_vcf_index = input_vcf + ".tbi"
  String base_name = basename(basename(basename(input_vcf, ".vcf.gz"), ".vcf.bgz"), ".bcf")
  String output_vcf = base_name + ".QC_ANNOTATED.vcf.gz"
  String output_tbi = base_name + ".QC_ANNOTATED.vcf.gz.tbi"
  Int disk_size = ceil(size(input_vcf, 'GB') * 3) + 20

  command <<<
  set -euo
  
  input_file="~{input_vcf}"
  touch ~{input_vcf_index}
  # Use nproc-1 to leave buffer for system overhead
  CHUNKS=$(( $(nproc) - 1 ))
  if [ $CHUNKS -lt 1 ]; then CHUNKS=1; fi
  
  echo "=== Parallel Filter by Chromosome ==="
  echo "Input: $input_file"
  echo "Output: ~{output_vcf}"
  echo "Genotype filter: ~{genotype_filter}"
  echo "Variant filter: ~{variant_filter}"
  echo "CPU cores: $CHUNKS"
  echo ""
  
  # Get list of chromosomes from the VCF index
  echo "Extracting chromosome list..."
  mapfile -t chromosomes < <(bcftools index -s "$input_file" | awk '{print $1}')
  num_chroms=${#chromosomes[@]}
  
  echo "Found $num_chroms chromosomes: ${chromosomes[*]}"
  
  if [[ $num_chroms -eq 0 ]]; then
      echo "Error: No chromosomes found in VCF"
      exit 1
  fi

  # Create chromosome list file for parallel processing
  echo "Processing $num_chroms chromosomes in parallel (max $CHUNKS jobs)..."
  echo ""
  
  # Create a file with complete bcftools commands for each chromosome
  rm -f commands.txt
  for chrom in "${chromosomes[@]}"; do
    safe_chrom=$(echo "$chrom" | sed 's/[*:\/]/_/g')
    # Escape curly braces for eval
    cat >> commands.txt << EOF
echo "Processing: $chrom" && bcftools view -r "\\{$chrom\\}" "$input_file" -Ou | bcftools +setGT -Ou -- -t q -n . -i '~{genotype_filter}' | bcftools +fill-tags -Ou -- -t AC | bcftools view -i '~{variant_filter}' -Ou | bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' -Oz -o "chunk_${safe_chrom}.vcf.gz" && echo "  ✓ Done: $chrom"
EOF
  done
  
  echo "Sample commands (first 5):"
  head -5 commands.txt
  echo ""
  echo "Sample HLA commands:"
  grep "HLA" commands.txt | head -3 || echo "No HLA chromosomes found"
  echo ""
  
  # Execute all commands in parallel using bash background jobs
  while read -r cmd; do
    while [ $(jobs -r | wc -l) -ge $CHUNKS ]; do sleep 0.1; done
    eval "$cmd" &
  done < commands.txt
  wait
  
  echo ""
  echo "Concatenating chromosomes in original order..."
  # Build ordered list of chunk files based on chromosome order
  rm -f chunk_list.txt
  for chrom in "${chromosomes[@]}"; do
    safe_chrom=$(echo "$chrom" | sed 's/[*:\/]/_/g')
    echo "chunk_${safe_chrom}.vcf.gz" >> chunk_list.txt
  done
  
  # Concatenate in chromosome order
  bcftools concat -n -Oz -o ~{output_vcf} -f chunk_list.txt
  
  echo "Indexing final output..."
  tabix -p vcf ~{output_vcf}

  # Cleanup
  echo "Cleaning up temporary files..."
  rm -f chunk_*.vcf.gz commands.txt chunk_list.txt
  echo "=== Complete! ==="
  >>>

  output {
    File filtered_vcf = output_vcf
    File filtered_vcf_tbi = output_tbi
  }

  runtime {
    disks: "local-disk ~{disk_size} HDD"
    cpu: cpu_count
    preemptible: 1
  }
}

task ValidateFiltering {
  input {
    File original_vcf
    File filtered_vcf
    File filtered_vcf_tbi
    String genotype_filter
  }

  Int disk_size = ceil(size(original_vcf, 'GB') + size(filtered_vcf, 'GB')) + 10

  command <<<
    set -euo

    echo "=== Validating VCF Filtering ===" > report.txt
    echo "Original: ~{original_vcf}" >> report.txt
    echo "Filtered: ~{filtered_vcf}" >> report.txt
    echo "" >> report.txt

    # Touch index
    touch ~{filtered_vcf_tbi}

    # Test 1: Check genotypes were set to missing
    echo "Test 1: Checking genotypes were set to missing (~{genotype_filter})" >> report.txt
    echo "--------------------------------------------------------------" >> report.txt

    # Find genotypes that meet criteria and check if variant still exists in filtered file
    found=0
    bcftools query -i 'GT!="mis"' -f '%CHROM\t%POS[\t%SAMPLE\t%GT\t%DP\t%GQ]\n' "~{original_vcf}" 2>/dev/null | \
      awk '($4 < 10 || $5 < 20)' | \
      while read chrom pos sample gt dp gq; do
        # Check if this variant exists in filtered file
        if bcftools view -H -r "$chrom:$pos" "~{filtered_vcf}" 2>/dev/null | grep -q .; then
          after_gt=$(bcftools query -s "$sample" -r "$chrom:$pos" -f '[\t%GT]\n' "~{filtered_vcf}" 2>/dev/null | tr -d '\t' | tr -d ' ')
          
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

    first_var=$(bcftools view -H "~{filtered_vcf}" | head -1 | awk '{print $1":"$2}')

    if [[ -n "$first_var" ]]; then
      ac_value=$(bcftools query -r "$first_var" -f '%AC\n' "~{filtered_vcf}" 2>/dev/null)
      
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

    sample_id=$(bcftools view -H "~{filtered_vcf}" 2>/dev/null | head -1 | awk '{print $3}')
    if [[ "$sample_id" =~ ^[^_]+_[0-9]+_.+_.+$ ]]; then
      echo "✓ PASS: IDs formatted as CHROM_POS_REF_ALT (example: $sample_id)" >> report.txt
    else
      echo "✗ FAIL: ID not in expected format (got: $sample_id)" >> report.txt
    fi

    echo "" >> report.txt

    # Test 4: Check chromosome order preservation and filtering statistics
    echo "Test 4: Chromosome order and filtering statistics" >> report.txt
    echo "---------------------------------------------------" >> report.txt
    
    # Extract unique chromosomes from both files with counts
    bcftools query -f '%CHROM\n' "~{original_vcf}" | sort | uniq -c | awk '{print $2"\t"$1}' > original_chrom_counts.txt
    bcftools query -f '%CHROM\n' "~{filtered_vcf}" | sort | uniq -c | awk '{print $2"\t"$1}' > filtered_chrom_counts.txt
    
    # Get ordered chromosome list
    bcftools query -f '%CHROM\n' "~{original_vcf}" | uniq > original_chroms.txt
    bcftools query -f '%CHROM\n' "~{filtered_vcf}" | uniq > filtered_chroms.txt
    
    echo "Original chromosomes: $(wc -l < original_chroms.txt)" >> report.txt
    echo "Filtered chromosomes: $(wc -l < filtered_chroms.txt)" >> report.txt
    echo "" >> report.txt
    
    # Per-chromosome filtering statistics
    echo "Per-chromosome variant counts and filtering:" >> report.txt
    printf "%-30s %10s %10s %10s\n" "CHROMOSOME" "ORIGINAL" "FILTERED" "%_DROPPED" >> report.txt
    printf "%-30s %10s %10s %10s\n" "----------" "--------" "--------" "---------" >> report.txt
    
    while read chrom; do
      orig_count=$(grep -w "^$chrom" original_chrom_counts.txt | awk '{print $2}' || echo "0")
      filt_count=$(grep -w "^$chrom" filtered_chrom_counts.txt | awk '{print $2}' || echo "0")
      
      if [[ "$orig_count" -gt 0 ]]; then
        pct_dropped=$(awk "BEGIN {printf \"%.1f\", (($orig_count - $filt_count) / $orig_count) * 100}")
      else
        pct_dropped="0.0"
      fi
      
      printf "%-30s %10s %10s %9s%%\n" "$chrom" "$orig_count" "$filt_count" "$pct_dropped" >> report.txt
    done < original_chroms.txt
    
    echo "" >> report.txt
    
    # Check if filtered chromosomes are in same order
    if diff original_chroms.txt filtered_chroms.txt > /dev/null 2>&1; then
      echo "✓ PASS: All chromosomes preserved in original order" >> report.txt
    elif comm -12 original_chroms.txt filtered_chroms.txt | diff - filtered_chroms.txt > /dev/null 2>&1; then
      echo "✓ PASS: Chromosome order preserved (some chromosomes filtered out)" >> report.txt
    else
      echo "✗ FAIL: Chromosome order differs from input!" >> report.txt
    fi
    
    rm -f original_chroms.txt filtered_chroms.txt original_chrom_counts.txt filtered_chrom_counts.txt

    echo "" >> report.txt
    echo "=== Validation Complete ===" >> report.txt

    cat report.txt
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
