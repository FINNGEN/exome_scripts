version 1.0

workflow adpkd_filter {
  input {
    Array[File] vcf_files
    String genotype_filter
    String variant_filter
    Int cpu_count = 8
  }

  scatter (vcf in vcf_files) {
    File vcf_index = vcf + ".tbi"

    call FilterByChromosome {
      input:
        input_vcf = vcf,
        input_vcf_index = vcf_index,
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
    File input_vcf_index
    String genotype_filter
    String variant_filter
    Int cpu_count
  }

  Int disk_size = ceil(size(input_vcf, 'GB') * 3) + 20

  command <<<
    set -euo

    input_file="~{input_vcf}"
    CHUNKS=~{cpu_count}
    
    # Create output filename with QC_ANNOTATED suffix
    basename=$(basename "$input_file" .vcf.gz)
    basename=$(basename "$basename" .vcf.bgz)
    basename=$(basename "$basename" .bcf)
    output_file="${basename}.QC_ANNOTATED.vcf.gz"
    
    echo "=== ADPKD Parallel Filter by Chromosome ==="
    echo "Input: $input_file"
    echo "Output: $output_file"
    echo "Genotype filter: ~{genotype_filter}"
    echo "Variant filter: ~{variant_filter}"
    echo "CPU cores: $CHUNKS"
    echo ""

    # Touch index to ensure it's available
    touch ~{input_vcf_index}

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
    rm -f chrom_list.txt
    for chrom in "${chromosomes[@]}"; do
      echo "$chrom" >> chrom_list.txt
    done

    echo "Processing $num_chroms chromosomes in parallel (max $CHUNKS jobs)..."

    # Process each chromosome in parallel
    cat chrom_list.txt | parallel -j "$CHUNKS" '
      chrom={}
      output="chunk_{}.vcf.gz"
      echo "Processing chromosome: $chrom"
      
      bcftools view -r "$chrom" "'"$input_file"'" -Ou | \
        bcftools +setGT -Ou -- -t q -n . -i "~{genotype_filter}" | \
        bcftools +fill-tags -Ou -- -t AC | \
        bcftools view -i "~{variant_filter}" -Ou | \
        bcftools annotate --set-id +"%CHROM\_%POS\_%REF\_%ALT" -Oz -o "$output"
      
      echo "Completed chromosome: $chrom"
    '

    echo "Indexing chromosome chunks..."
    ls chunk_*.vcf.gz | parallel -j "$CHUNKS" 'tabix -p vcf {}'

    echo "Concatenating chromosomes..."
    bcftools concat -n -Oz -o filtered.vcf.gz chunk_*.vcf.gz

    echo "Indexing final output..."
    tabix -p vcf filtered.vcf.gz

    # Cleanup
    echo "Cleaning up temporary files..."
    rm -f chunk_*.vcf.gz chunk_*.vcf.gz.tbi chrom_list.txt

    echo "=== Complete! ==="
  >>>

  output {
    File filtered_vcf = "filtered.vcf.gz"
    File filtered_vcf_tbi = "filtered.vcf.gz.tbi"
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
