version 1.0

workflow single_file_qc {
  input {
    Array[File] vcf_files
    String genotype_filter
    String variant_filter
    Int cpu_count
    Int? test_sample_count
    File norm_fasta
  }

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
        input_vcf = vcf_to_filter
    }

    call FilterByChromosome {
      input:
        input_vcf = vcf_to_filter,
        genotype_filter = genotype_filter,
        variant_filter = variant_filter,
        cpu_count = cpu_count,
        norm_fasta = norm_fasta
    }

    call ComputeStats as FilteredStats {
      input:
        input_vcf = FilterByChromosome.filtered_vcf
    }

    call ValidateFiltering {
      input:
        original_vcf = vcf_to_filter,
        original_chrom_counts = OriginalStats.chrom_counts,
        filtered_vcf = FilterByChromosome.filtered_vcf,
        filtered_vcf_tbi = FilterByChromosome.filtered_vcf_tbi,
        filtered_chrom_counts = FilteredStats.chrom_counts,
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
    File norm_fasta
  }

  File input_vcf_index = input_vcf + ".tbi"
  File norm_fasta_fai = norm_fasta + ".fai"
  String base_name = basename(basename(basename(input_vcf, ".vcf.gz"), ".vcf.bgz"), ".bcf")
  String output_vcf = base_name + ".QC_ANNOTATED.vcf.gz"
  String output_tbi = base_name + ".QC_ANNOTATED.vcf.gz.tbi"
  Int disk_size = ceil(size(input_vcf, 'GB') * 3) + 20
  Int vcf_size_gb = ceil(size(input_vcf, 'GB'))
  Int memory_gb = if vcf_size_gb < 8 then 8 else vcf_size_gb

  command <<<
  set -euo

  input_file="~{input_vcf}"
  touch ~{input_vcf_index}
  touch ~{norm_fasta_fai}
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
    cat >> commands.txt << EOF
echo "Processing: $chrom" && bcftools view -r "{$chrom}" "$input_file" -Ou | bcftools norm -f '~{norm_fasta}' -c x -Ou | bcftools +setGT -Ou -- -t q -n . -i '~{genotype_filter}' | bcftools +fill-tags -Ou -- -t AC | bcftools view -i '~{variant_filter}' -Ou | bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' -Oz -o "chunk_${safe_chrom}.vcf.gz" && echo "  ✓ Done: $chrom"
EOF
  done
  
  echo "Sample commands (first 5):"
  head -5 commands.txt
  echo ""
  echo "Sample HLA commands:"
  grep "HLA" commands.txt | head -3 || echo "No HLA chromosomes found"
  echo ""
  
  # Execute all commands in parallel
  parallel -j "$CHUNKS" < commands.txt
  
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
    memory: "~{memory_gb}G"
    disks: "local-disk ~{disk_size} HDD"
    cpu: cpu_count
    preemptible: 1
  }
}

task ComputeStats {
  input {
    File input_vcf
  }

  Int disk_size = ceil(size(input_vcf, 'GB')) + 10

  command <<<
  set -euo

  echo "=== Computing chromosome statistics ==="

  # Extract unique chromosomes with counts
  bcftools query -f '%CHROM\n' "~{input_vcf}" | sort | uniq -c | awk '{print $2"\t"$1}' > chrom_counts.txt

  echo "Found $(wc -l < chrom_counts.txt) chromosomes"
  >>>

  output {
    File chrom_counts = "chrom_counts.txt"
  }

  runtime {
    memory: "4G"
    disks: "local-disk ~{disk_size} HDD"
    cpu: 1
    preemptible: 1
  }
}

task ValidateFiltering {
  input {
    File original_vcf
    File original_chrom_counts
    File filtered_vcf
    File filtered_vcf_tbi
    File filtered_chrom_counts
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
          echo "✓ PASS: $chrom:$pos GT: $gt → $after_gt" >> report.txt
        else
          echo "✗ FAIL: $chrom:$pos GT: $gt → $after_gt (not missing)" >> report.txt
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

  # Test 4: Filtering statistics by chromosome
  echo "Test 4: Filtering statistics by chromosome" >> report.txt
  echo "-------------------------------------------" >> report.txt

  # Use pre-computed stats for both files
  cp "~{original_chrom_counts}" original_chrom_counts.txt
  cp "~{filtered_chrom_counts}" filtered_chrom_counts.txt

  echo "Original chromosomes: $(wc -l < original_chrom_counts.txt)" >> report.txt
  echo "" >> report.txt

  # Per-chromosome filtering statistics (sorted by original count descending)
  echo "Per-chromosome variant counts (sorted by original count):" >> report.txt
  printf "%-30s %10s %10s %10s\n" "CHROMOSOME" "ORIGINAL" "FILTERED" "%_DROPPED" >> report.txt
  printf "%-30s %10s %10s %10s\n" "----------" "--------" "--------" "---------" >> report.txt

  # Sort chromosomes by original count (descending)
  sort -t$'\t' -k2 -nr original_chrom_counts.txt | while IFS=$'\t' read chrom orig_count; do
    filt_count=$(grep -Fw "$chrom" filtered_chrom_counts.txt | awk '{print $2}')
    [[ -z "$filt_count" ]] && filt_count=0
    
    if [[ "$orig_count" -gt 0 ]]; then
      pct_dropped=$(awk "BEGIN {printf \"%.1f\", (($orig_count - $filt_count) / $orig_count) * 100}")
    else
      pct_dropped="0.0"
    fi
    
    printf "%-30s %10s %10s %9s%%\n" "$chrom" "$orig_count" "$filt_count" "$pct_dropped" >> report.txt
  done

  echo "" >> report.txt

  echo "" >> report.txt
  echo "=== Validation Complete ===" >> report.txt

  cat report.txt
  >>>

  output {
    File report = "report.txt"
  }

  runtime {
    disks: "local-disk ~{disk_size} HDD"
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

