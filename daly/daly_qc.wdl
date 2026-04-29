version 1.0

workflow daly_qc {
  input {
    File vcf_list
    Boolean test_mode
  }
  Array[File] vcf_files = read_lines(vcf_list)
  scatter (vcf in vcf_files) {
    call AnnotateHeaders {
      input:
        input_vcf = vcf,
        test_mode = test_mode
    }
    
    call ParallelFilter {
      input:
        input_vcf = AnnotateHeaders.annotated_vcf,
        input_vcf_tbi = AnnotateHeaders.annotated_vcf_tbi
    }
  }

  call SortAndMerge {
    input:
      vcf_files = ParallelFilter.filtered_vcf,
      vcf_tbi_files = ParallelFilter.filtered_vcf_tbi
  }

  output {
    Array[File] filtered_vcfs = ParallelFilter.filtered_vcf
    Array[File] filtered_vcf_tbis = ParallelFilter.filtered_vcf_tbi
    File merged_vcf = SortAndMerge.merged_vcf
    File merged_vcf_tbi = SortAndMerge.merged_vcf_tbi
  }
}

task AnnotateHeaders {
  input {
    File input_vcf
    Boolean test_mode
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
  
  # Select (possibly temporary) VCF to operate on
  vcf_for_annotation="~{input_vcf}"
  if [[ "~{test_mode}" == "true" ]]; then
      # Use only first 10k variants
      tmp_test_file=$(mktemp --suffix=.vcf.gz)
      (bcftools view -h "~{input_vcf}"; bcftools view -H "~{input_vcf}" | head -n 10000) | bgzip > "$tmp_test_file"
      tabix -p vcf "$tmp_test_file"
      vcf_for_annotation="$tmp_test_file"
  fi

  # Extract FILTER IDs in header (populate array)
  mapfile -t header_filters < <(bcftools view -h "$vcf_for_annotation" | awk -F'[=,]' '/^##FILTER=/{print $3}')
  
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
      bcftools annotate --header-lines "$missing_header" "$vcf_for_annotation" -Oz -o "$output_vcf"
  else
      cp "$vcf_for_annotation"  "$output_vcf"
  fi
  tabix -p vcf "$output_vcf"
  
  >>>

  output {
    File annotated_vcf = "annotated.${basename(input_vcf)}"
    File annotated_vcf_tbi = "annotated.${basename(input_vcf)}.tbi"
  }

  runtime {
    memory: "4G"
    disks: "local-disk ~{disk_size} HDD"
    cpu: 1
    preemptible: 1
  }
}

task ParallelFilter {
  input {
    File input_vcf
    File input_vcf_tbi
  }

  Int disk_size = ceil(size(input_vcf,'GB')*3) + 20
  Int cpu_count = 8
  
  command <<<
  set -euo pipefail
  
  input_file="~{input_vcf}"
  CHUNKS=~{cpu_count}

  echo "Creating $CHUNKS region files..."
  # Get chromosome, first position, and last position
  read chrom first_pos last_pos < <(bcftools query -f'%CHROM\t%POS\n' "$input_file" | awk 'NR==1{chrom=$1; first=$2} {last=$2} END{print chrom, first, last}')

  # Split interval equally using Python linspace
  python3 -c "import numpy as np; [open(f'region_chunk_{i:02d}','w').write(f'$chrom\t{int(s)}\t{int(e)}\n') for i,(s,e) in enumerate(zip(np.linspace($first_pos,$last_pos,$CHUNKS+1)[:-1], np.linspace($first_pos,$last_pos,$CHUNKS+1)[1:]))]"

  echo "Processing chunks in parallel..."
  ls region_chunk_* | parallel -j $CHUNKS 'bcftools view '"$input_file"' -R {} -Ou | bcftools view -e '"'"'FILTER~"NO_HQ_GENOTYPES"'"'"' -Ou | bcftools annotate --set-id +'"'"'%CHROM\_%POS\_%REF\_%ALT'"'"' -Oz -o chunk{}.vcf.gz && echo {}'

  echo "Concatenating chunks..."
  bcftools concat -n -f <(ls chunk*.vcf.gz) -Oz -o filtered.vcf.gz && rm chunk*.vcf.gz region_chunk_*

  echo "Indexing output..."
  tabix -p vcf filtered.vcf.gz
  
  echo "Done! Output file: filtered.vcf.gz"
  >>>

  output {
    File filtered_vcf = "filtered.vcf.gz"
    File filtered_vcf_tbi = "filtered.vcf.gz.tbi"
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
  }

  Int disk_size = ceil(size(vcf_files,'GB')*2) + 50
  
  command <<<
  set -euo pipefail
  
  # Create file list
  cat ~{write_lines(vcf_files)} > unsorted_vcf_list.txt
  
  # Sort by chromosome (extract chr from filename and sort naturally)
  sort -V unsorted_vcf_list.txt > sorted_vcf_list.txt
  
  echo "Sorted VCF files:"
  cat sorted_vcf_list.txt
  
  echo "Concatenating sorted chromosomes..."
  bcftools concat -f sorted_vcf_list.txt -Oz -o merged.vcf.gz
  
  echo "Indexing merged VCF..."
  tabix -p vcf merged.vcf.gz
  
  echo "Done! Output file: merged.vcf.gz"
  >>>

  output {
    File merged_vcf = "merged.vcf.gz"
    File merged_vcf_tbi = "merged.vcf.gz.tbi"
  }

  runtime {
    memory: "8G"
    disks: "local-disk ~{disk_size} HDD"
    cpu: 2
    preemptible: 1
  }
}
