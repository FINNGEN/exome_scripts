version 1.0

workflow single_file_qc {
  input {
    Array[File] vcf_files
    String genotype_filter
    String variant_filter
    Int cpu_count
    Int? test_sample_count
    File norm_fasta
    File denials                  # sample IDs to remove, one per line (see FilterByChromosome)
    File aliases                  # tab-delimited alias groups (one group per line, same members
                                   # share a group) — used to expand `denials` to every alias of
                                   # each denied ID before exclusion, since the ID actually present
                                   # in a given VCF's header may be an alias rather than the ID
                                   # recorded in the denials list (see ExpandDenials)
  }

  # Expand denials to cover alias variants before any per-chromosome FilterByChromosome call uses it.
  # Dataset-independent, so this runs once for the whole workflow rather than once per VCF.
  call ExpandDenials {
    input:
      denials = denials,
      aliases = aliases
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
        norm_fasta = norm_fasta,
        denials = ExpandDenials.expanded_denials
    }

    call ComputeStats as FilteredStats {
      input:
        input_vcf = FilterByChromosome.filtered_vcf
    }

    call ValidateFiltering {
      input:
        filtered_vcf_name     = basename(FilterByChromosome.filtered_vcf),
        original_sample_vcf   = OriginalStats.sample_vcf,
        original_chrom_counts = OriginalStats.chrom_counts,
        filtered_sample_vcf   = FilteredStats.sample_vcf,
        filtered_chrom_counts = FilteredStats.chrom_counts,
        genotype_filter       = genotype_filter
    }
  }

  output {
    Array[File] filtered_vcfs = FilterByChromosome.filtered_vcf
    Array[File] filtered_vcf_tbis = FilterByChromosome.filtered_vcf_tbi
    Array[File] validation_reports = ValidateFiltering.report
  }
}

# ── ExpandDenials ─────────────────────────────────────────────────────────────
# Expands a raw denials list (one ID per line) to include every alias of each
# denied ID, per the alias groups in `aliases` (same tab-delimited group-per-
# line format as resolve_mapping.py's DEFAULT_ALIASES / load_aliases). Needed
# because the ID recorded in the denials list and the ID actually present in a
# given VCF's header can be different aliases of the same participant —
# bcftools view -S ^denials only matches literal strings, so
# FilterByChromosome would otherwise fail to exclude (or hard-error on) a
# denied sample whose header ID is an alias rather than the one on the list.

task ExpandDenials {
  input {
    File denials    # plain list, one ID per line
    File aliases    # tab-delimited alias groups, one group per line
  }

  command <<<
  set -euo pipefail
  python3 << 'PY'
denied = {line.strip() for line in open("~{denials}") if line.strip()}
expanded = set(denied)
with open("~{aliases}") as fh:
    for line in fh:
        ids = [x.strip() for x in line.strip().split("\t") if x.strip()]
        if len(ids) < 2:
            continue
        if denied & set(ids):
            expanded.update(ids)
with open("expanded_denials.txt", "w") as out:
    for id_ in sorted(expanded):
        out.write(id_ + "\n")
print(f"Expanded {len(denied)} denied IDs to {len(expanded)} IDs (incl. aliases)")
PY
  >>>

  output {
    File expanded_denials = "expanded_denials.txt"
  }

  runtime {
    memory: "2G"
    disks: "local-disk 10 HDD"
    cpu: 1
    preemptible: 1
  }
}

task FilterByChromosome {
  input {
    File input_vcf
    String genotype_filter
    String variant_filter
    Int cpu_count
    File norm_fasta
    File denials    # sample IDs to remove, one per line
  }

  File input_vcf_index = input_vcf + ".tbi"
  File norm_fasta_fai = norm_fasta + ".fai"
  String base_name = basename(basename(basename(input_vcf, ".vcf.gz"), ".vcf.bgz"), ".bcf")
  String output_vcf = base_name + ".QC_ANNOTATED.vcf.gz"
  String output_tbi = base_name + ".QC_ANNOTATED.vcf.gz.tbi"
  Int disk_size = ceil(size(input_vcf, 'GB') * 3) + 20
  Int memory_gb = cpu_count * 2 + 4

  command <<<
  set -euo
  INPUT_VCF="~{input_vcf}"
  touch "~{input_vcf_index}"
  touch "~{norm_fasta_fai}"
  NORM_FASTA="~{norm_fasta}"
  GENOTYPE_FILTER='~{genotype_filter}'
  VARIANT_FILTER='~{variant_filter}'
  OUTPUT_VCF="~{output_vcf}"
  DENIALS="~{denials}"
  chromosomes=()
  CHUNKS=$(( $(nproc) - 1 ))
  if [[ $CHUNKS -lt 1 ]]; then CHUNKS=1; fi

  echo "=== Parallel Filter by Chromosome ==="
  echo "Input:           $INPUT_VCF"
  echo "Output:          $OUTPUT_VCF"
  echo "Genotype filter: $GENOTYPE_FILTER"
  echo "Variant filter:  $VARIANT_FILTER"
  echo "CPU cores:       $CHUNKS"
  echo "Removing denied samples listed in: $DENIALS"
  echo ""

  if [[ ${#chromosomes[@]} -eq 0 ]]; then
    mapfile -t chromosomes < <(bcftools index -s "$INPUT_VCF" | awk '{print $1}')
  fi
  num_chroms=${#chromosomes[@]}
  echo "Found $num_chroms chromosomes: ${chromosomes[*]}"

  if [[ $num_chroms -eq 0 ]]; then
    echo "Error: No chromosomes found in VCF"
    exit 1
  fi

  echo "Processing $num_chroms chromosomes in parallel (max $CHUNKS jobs)..."
  echo ""

  # Detect chromosome naming: if VCF uses non-chr names, rename to chr prefix in output
  FIRST_CHROM=$(bcftools index -s "$INPUT_VCF" | awk 'NR==1{print $1}')
  RENAME_TO_CHR=""
  if [[ "$FIRST_CHROM" != chr* ]]; then
    echo "VCF uses non-chr chromosome names — renaming to chr prefix in output"
    RENAME_TO_CHR=$(mktemp)
    for i in $(seq 1 22) X Y MT M; do
      echo "$i chr$i" >> "$RENAME_TO_CHR"
    done
  fi
  RENAME_STEP="cat"
  [[ -n "$RENAME_TO_CHR" ]] && RENAME_STEP="bcftools annotate --rename-chrs ${RENAME_TO_CHR} -Ou"

  # Generate per-chromosome scripts and run in parallel
  # Variables expand at generation time — no quoting or function-export issues
  SCRIPT_DIR=$(mktemp -d)
  for chrom in "${chromosomes[@]}"; do
    safe_chrom=$(echo "$chrom" | sed 's/[*:\/]/_/g')
    output="chunk_${safe_chrom}.vcf.gz"
    cat > "${SCRIPT_DIR}/run_${safe_chrom}.sh" << SCRIPT
#!/bin/bash
echo "Processing chromosome: ${chrom}"
printf '${chrom}\t0\t9999999999\n' > "${output}.region.bed"
bcftools view -S ^${DENIALS} --force-samples -R "${output}.region.bed" "${INPUT_VCF}" | \\
    tr -d '\0' | \\
    ${RENAME_STEP} | \\
    bcftools norm -f '${NORM_FASTA}' -m -any -c x -Ou | \\
    bcftools +setGT -Ou -- -t q -n . -i '${GENOTYPE_FILTER}' | \\
    bcftools +fill-tags -Ou -- -t AC | \\
    bcftools view -i '${VARIANT_FILTER}' -Ou | \\
    bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' -Oz -o "${output}"
rm -f "${output}.region.bed"
echo "Completed chromosome: ${chrom}"
SCRIPT
  done

  ls "${SCRIPT_DIR}"/run_*.sh | parallel -j "$CHUNKS" 'bash {}'
  rm -rf "$SCRIPT_DIR"

  echo ""
  echo "Concatenating chromosomes in original order..."
  rm -f chunk_list.txt
  for chrom in "${chromosomes[@]}"; do
    safe_chrom=$(echo "$chrom" | sed 's/[*:\/]/_/g')
    echo "chunk_${safe_chrom}.vcf.gz" >> chunk_list.txt
  done

  bcftools concat -n -Oz -o "$OUTPUT_VCF" -f chunk_list.txt

  echo "Indexing final output..."
  tabix -p vcf "$OUTPUT_VCF"

  # Cleanup
  rm -f chunk_*.vcf.gz chunk_list.txt
  [[ -n "$RENAME_TO_CHR" ]] && rm -f "$RENAME_TO_CHR"
  echo "=== Complete! ==="
  >>>

  output {
    File filtered_vcf = output_vcf
    File filtered_vcf_tbi = output_tbi
  }

  runtime {
    memory: "~{memory_gb} GB"
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
  touch "~{input_vcf_tbi}"

  echo "=== Computing chromosome statistics ==="
  bcftools index -s "~{input_vcf}" | awk '{print $1"\t"$3}' > chrom_counts.txt
  echo "Found $(wc -l < chrom_counts.txt) chromosomes"

  echo "Creating sample VCF for validation (100 variants)..."
  bcftools view -h "~{input_vcf}" | bgzip -c > sample.vcf.gz
  bcftools view -H "~{input_vcf}" | head -n 100 | bgzip -c >> sample.vcf.gz
  tabix -p vcf sample.vcf.gz
  echo "Sample VCF created: $(bcftools view -H sample.vcf.gz | wc -l) variants"
  >>>

  output {
    File chrom_counts   = "chrom_counts.txt"
    File sample_vcf     = "sample.vcf.gz"
    File sample_vcf_tbi = "sample.vcf.gz.tbi"
  }

  runtime {
    memory: "4 GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: 1
  }
}

task ValidateFiltering {
  input {
    String filtered_vcf_name
    File original_sample_vcf
    File original_chrom_counts
    File filtered_sample_vcf
    File filtered_chrom_counts
    String genotype_filter
  }

  File   original_sample_vcf_tbi = original_sample_vcf + ".tbi"
  File   filtered_sample_vcf_tbi = filtered_sample_vcf + ".tbi"
  String report_name             = sub(filtered_vcf_name, "\\.vcf\\.gz$", ".report.txt")
  Int    disk_size               = ceil(size(original_sample_vcf, 'GB') + size(filtered_sample_vcf, 'GB')) + 10

  command <<<
  set -euo

  echo "=== Validating VCF Filtering ===" > ~{report_name}
  echo "Original: ~{original_sample_vcf}" >> ~{report_name}
  echo "Filtered: ~{filtered_sample_vcf}" >> ~{report_name}
  echo "" >> ~{report_name}

  # Touch index
  touch ~{filtered_sample_vcf_tbi}

  # Test 1: Check genotypes were set to missing
  echo "Test 1: Checking genotypes were set to missing (~{genotype_filter})" >> ~{report_name}
  echo "--------------------------------------------------------------" >> ~{report_name}

  # Find genotypes that meet criteria and check if variant still exists in filtered file
  found=0
  bcftools query -i 'GT!="mis"' -f '%CHROM\t%POS[\t%SAMPLE\t%GT\t%DP\t%GQ]\n' "~{original_sample_vcf}" 2>/dev/null | \
    awk '($4 < 10 || $5 < 20)' | \
    while read chrom pos sample gt dp gq; do
      # Check if this variant exists in filtered file
      if bcftools view -H -r "$chrom:$pos" "~{filtered_sample_vcf}" 2>/dev/null | grep -q .; then
        after_gt=$(bcftools query -s "$sample" -r "$chrom:$pos" -f '[\t%GT]\n' "~{filtered_sample_vcf}" 2>/dev/null | tr -d '\t' | tr -d ' ')
        
        if [[ "$after_gt" == "./." ]] || [[ "$after_gt" == "." ]]; then
          echo "✓ PASS: $chrom:$pos GT: $gt → $after_gt" >> ~{report_name}
        else
          echo "✗ FAIL: $chrom:$pos GT: $gt → $after_gt (not missing)" >> ~{report_name}
        fi
        found=1
        break
      fi
    done

  if [[ $found -eq 0 ]]; then
    echo "⚠ NOTE: All low-quality genotypes resulted in variants being removed (AC=0)" >> ~{report_name}
  fi

  echo "" >> ~{report_name}

  # Test 2: Check AC field exists
  echo "Test 2: Checking AC field is present in filtered VCF" >> ~{report_name}
  echo "-----------------------------------------------------" >> ~{report_name}

  first_var=$(bcftools view -H "~{filtered_sample_vcf}" | head -1 | awk '{print $1":"$2}')

  if [[ -n "$first_var" ]]; then
    ac_value=$(bcftools query -r "$first_var" -f '%AC\n' "~{filtered_sample_vcf}" 2>/dev/null)
    
    if [[ -n "$ac_value" ]]; then
      echo "✓ PASS: AC field present (example: $first_var AC=$ac_value)" >> ~{report_name}
    else
      echo "✗ FAIL: AC field missing in filtered VCF" >> ~{report_name}
    fi
  else
    echo "⚠ NOTE: No variants in filtered file" >> ~{report_name}
  fi

  echo "" >> ~{report_name}

  # Test 3: Check variant IDs
  echo "Test 3: Checking variant IDs (CHROM_POS_REF_ALT format)" >> ~{report_name}
  echo "--------------------------------------------------------" >> ~{report_name}

  sample_id=$(bcftools view -H "~{filtered_sample_vcf}" 2>/dev/null | head -1 | awk '{print $3}')
  if [[ "$sample_id" =~ ^[^_]+_[0-9]+_.+_.+$ ]]; then
    echo "✓ PASS: IDs formatted as CHROM_POS_REF_ALT (example: $sample_id)" >> ~{report_name}
  else
    echo "✗ FAIL: ID not in expected format (got: $sample_id)" >> ~{report_name}
  fi

  echo "" >> ~{report_name}

  # Test 4: Filtering statistics by chromosome
  echo "Test 4: Filtering statistics by chromosome" >> ~{report_name}
  echo "-------------------------------------------" >> ~{report_name}

  # Use pre-computed stats for both files
  cp "~{original_chrom_counts}" original_chrom_counts.txt
  cp "~{filtered_chrom_counts}" filtered_chrom_counts.txt

  echo "Original chromosomes: $(wc -l < original_chrom_counts.txt)" >> ~{report_name}
  echo "" >> ~{report_name}

  # Per-chromosome filtering statistics (sorted by original count descending)
  echo "Per-chromosome variant counts (sorted by original count):" >> ~{report_name}
  printf "%-30s %10s %10s %10s\n" "CHROMOSOME" "ORIGINAL" "FILTERED" "%_DROPPED" >> ~{report_name}
  printf "%-30s %10s %10s %10s\n" "----------" "--------" "--------" "---------" >> ~{report_name}

  # Sort chromosomes by original count (descending)
  # Compare chrom names with any "chr" prefix stripped on both sides: FilterByChromosome
  # renames to chr-prefixed output when the input isn't already chr-prefixed, so a literal
  # name join here would otherwise silently fail (e.g. "1" vs "chr1") and default every
  # filtered count to 0. Match on column 1 only (awk $1==c), not grep -Fw against the whole
  # line — with alt-contig chromosomes, a short normalized name like "1" can otherwise
  # collide with an unrelated line's *count* column when that count happens to equal "1".
  sed 's/^chr//' filtered_chrom_counts.txt > filtered_chrom_counts_norm.txt

  sort -t$'\t' -k2 -nr original_chrom_counts.txt | while IFS=$'\t' read chrom orig_count; do
    chrom_norm="${chrom#chr}"
    filt_count=$(awk -F'\t' -v c="$chrom_norm" '$1==c {print $2; exit}' filtered_chrom_counts_norm.txt)
    [[ -z "$filt_count" ]] && filt_count=0
    
    if [[ "$orig_count" -gt 0 ]]; then
      pct_dropped=$(awk "BEGIN {printf \"%.1f\", (($orig_count - $filt_count) / $orig_count) * 100}")
    else
      pct_dropped="0.0"
    fi
    
    printf "%-30s %10s %10s %9s%%\n" "$chrom" "$orig_count" "$filt_count" "$pct_dropped" >> ~{report_name}
  done

  # Totals across all chromosomes. Computed independently of the loop above (which runs in a
  # subshell via the pipe from `sort`, so its per-iteration variables don't survive past `done`).
  orig_total=$(awk -F'\t' '{s+=$2} END{print s+0}' original_chrom_counts.txt)
  filt_total=$(awk -F'\t' '{s+=$2} END{print s+0}' filtered_chrom_counts_norm.txt)
  if [[ "$orig_total" -gt 0 ]]; then
    total_pct_dropped=$(awk "BEGIN {printf \"%.1f\", (($orig_total - $filt_total) / $orig_total) * 100}")
  else
    total_pct_dropped="0.0"
  fi
  printf "%-30s %10s %10s %10s\n" "----------" "--------" "--------" "---------" >> ~{report_name}
  printf "%-30s %10s %10s %9s%%\n" "TOTAL" "$orig_total" "$filt_total" "$total_pct_dropped" >> ~{report_name}

  echo "" >> ~{report_name}

  # Test 5: Sample count (denial exclusion check)
  echo "Test 5: Sample count (denial exclusion check)" >> ~{report_name}
  echo "-----------------------------------------------" >> ~{report_name}

  orig_samples=$(bcftools query -l "~{original_sample_vcf}" | wc -l)
  filt_samples=$(bcftools query -l "~{filtered_sample_vcf}" | wc -l)
  samples_removed=$((orig_samples - filt_samples))

  echo "Original samples: $orig_samples" >> ~{report_name}
  echo "Filtered samples: $filt_samples" >> ~{report_name}
  echo "Samples removed:  $samples_removed" >> ~{report_name}

  echo "" >> ~{report_name}
  echo "=== Validation Complete ===" >> ~{report_name}

  cat ~{report_name}
  >>>

  output {
    File report = "~{report_name}"
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
    disks: "local-disk ~{disk_size} HDD"
  }
}

