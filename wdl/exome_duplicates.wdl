version 1.0

workflow exome_duplicates {
  input {
    Array[Array[String]] vcf_pairs
    File   plink_bed
    String plink_prefix
    File?  bim
  }

  File plink_bim = sub(plink_bed, "\\.bed$", ".bim")
  Array[File] plink_input_files = [plink_bed, plink_bim, sub(plink_bed, "\\.bed$", ".fam")]

  scatter (pair in vcf_pairs) {
    # pre-filter VCF to bim SNPs (parallel bcftools per chrom, streaming)
    call SubsetVCF {
      input:
        prefix    = pair[0],
        input_vcf = pair[1],
        bim       = select_first([bim, plink_bim])
    }
    # convert pre-filtered VCF to plink
    call ConvertToPlink as VcfToPlink {
      input:
        prefix      = pair[0],
        input_files = [SubsetVCF.filtered_vcf],
        bim         = select_first([bim, plink_bim])
    }
    # subset plink reference to the shared SNPs with exome data
    call ConvertToPlink as PlinkFilter {
      input:
        prefix      = plink_prefix,
        input_files = plink_input_files,
        bim         = VcfToPlink.plink_bim,
        memory_gb   = 32
    }

    # run KING to find duplicates between the exome dataset and the plink reference
    call RunKinship {
      input:
        vcf_plink = VcfToPlink.plink_data,
        ref_plink = PlinkFilter.plink_data,
        prefix    = pair[0] + "_king"
    }
  }

  output {
    Array[Array[File]] exome_plink = VcfToPlink.plink_data
    Array[File]        kinship_con = RunKinship.con_file
  }
}

task SubsetVCF {
  input {
    String prefix
    File   input_vcf
    File   bim
    Int    cpu       = 24
    Int    memory_gb = cpu
  }

  File   input_vcf_tbi = input_vcf + ".tbi"
  String output_vcf    = prefix + ".subset.vcf.gz"
  Int    disk_size     = ceil(size(input_vcf, 'GB') * 3)

  command <<<
  set -euo pipefail
  VCF="~{input_vcf}"
  touch "~{input_vcf_tbi}"
  BIM="~{bim}"
  OUTPUT_VCF="~{output_vcf}"
  CHUNKS=$(( $(nproc) - 1 ))
  if [[ $CHUNKS -lt 1 ]]; then CHUNKS=1; fi

  # chromosomes present in the bim, normalised to chr-prefix
  mkdir -p ./tmp
  mapfile -t chromosomes < <(awk '{print $1}' "$BIM" | sort -u | awk '{print (/^chr/ ? $0 : "chr"$0)}')

  N_BEFORE=$(bcftools index -s "$VCF" | awk '{sum+=$3} END {print sum}')
  echo "Variants before filter: $N_BEFORE"
  echo "Processing ${#chromosomes[@]} chromosomes in parallel ($CHUNKS jobs)..."
  echo ""

  # build per-chromosome CHROM\tPOS files; -T uses the tabix index for O(1) position lookup
  # rather than scanning the full VCF, so this is the main bottleneck reducer
  for chrom in "${chromosomes[@]}"; do
    safe_chrom=$(echo "$chrom" | sed 's/[*:\/]/_/g')
    awk -v c="$chrom" '$1 == c || "chr"$1 == c {print c"\t"$4}' "$BIM" > "./tmp/pos_${safe_chrom}.txt"
  done

  # generate and run one bcftools job per chromosome; ConvertToPlink does the final
  # ID-based --extract to handle any allele ambiguities after this positional pre-filter
  SCRIPT_DIR=$(mktemp -d)
  for chrom in "${chromosomes[@]}"; do
    safe_chrom=$(echo "$chrom" | sed 's/[*:\/]/_/g')
    cat > "${SCRIPT_DIR}/run_${safe_chrom}.sh" << SCRIPT
#!/bin/bash
bcftools view -r "${chrom}" -T "./tmp/pos_${safe_chrom}.txt" "${VCF}" -Oz -o "chunk_${safe_chrom}.vcf.gz"
echo "Done: ${chrom}"
SCRIPT
  done

  ls "${SCRIPT_DIR}"/run_*.sh | parallel -j "$CHUNKS" 'bash {}'

  # concat in original chromosome order (chunks are already sorted, -n skips re-sorting)
  for chrom in "${chromosomes[@]}"; do
    safe_chrom=$(echo "$chrom" | sed 's/[*:\/]/_/g')
    echo "chunk_${safe_chrom}.vcf.gz"
  done > chunk_list.txt

  bcftools concat -n -Oz -o "$OUTPUT_VCF" -f chunk_list.txt
  rm -f chunk_*.vcf.gz chunk_list.txt

  bcftools index -t "$OUTPUT_VCF"
  N_AFTER=$(bcftools index -s "$OUTPUT_VCF" | awk '{sum+=$3} END {print sum}')
  echo "Variants after filter: $N_AFTER"
  >>>

  output {
    File filtered_vcf     = output_vcf
    File filtered_vcf_tbi = output_vcf + ".tbi"
  }

  runtime {
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_size} HDD"
    cpu: cpu
  }
}

task RunKinship {
  input {
    Array[File] vcf_plink
    Array[File] ref_plink
    String      prefix
    Int         cpu       = 16
    Int         memory_gb = 64
  }

  Int disk_size = ceil(size(vcf_plink[0], 'GB') + size(ref_plink[0], 'GB')) * 2 + 10
  String docker = "eu.gcr.io/finngen-refinery-dev/exome_bioinf:king"

  command <<<
  set -euo pipefail

  VCF_BED="~{vcf_plink[0]}"
  REF_BED="~{ref_plink[0]}"
  VCF_PREFIX="${VCF_BED%.bed}"
  REF_PREFIX="${REF_BED%.bed}"
  OUTPUT="~{prefix}"

  echo "VCF dataset: $VCF_PREFIX ($(wc -l < "${VCF_PREFIX}.fam") samples, $(wc -l < "${VCF_PREFIX}.bim") SNPs)"
  echo "Ref dataset: $REF_PREFIX ($(wc -l < "${REF_PREFIX}.fam") samples, $(wc -l < "${REF_PREFIX}.bim") SNPs)"
  echo ""

  echo "Running KING --duplicate..."
  king -b "${VCF_BED}","${REF_BED}" --duplicate --prefix "$OUTPUT" --cpus $(nproc)

  if [[ ! -f "${OUTPUT}.con" ]]; then
    touch "${OUTPUT}.con"
  fi

  echo ""
  echo "Done."
  N_DUPS=$(tail -n +2 "${OUTPUT}.con" | wc -l)
  echo "  Duplicate pairs: ${N_DUPS}"
  if [[ $N_DUPS -gt 0 ]]; then
    echo ""
    head -5 "${OUTPUT}.con"
  fi
  >>>

  output {
    File con_file = "~{prefix}.con"
  }

  runtime {
    docker: docker
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_size} HDD"
    cpu: cpu
  }
}

task ConvertToPlink {
  input {
    String      prefix
    Array[File] input_files
    File        bim
    Int         memory_gb = 16
  }

  Int disk_size = ceil(size(input_files[0], 'GB') * 3) + 20

  command <<<
  set -euo
  THREADS=$(nproc)
  PREFIX="~{prefix}"
  BIM="~{bim}"
  INPUT_FILES=(~{sep=" " input_files})
  INPUT="${INPUT_FILES[0]}"

  cut -f2 "$BIM" > snp_list.txt
  SNP_LIST="snp_list.txt"
  echo "SNP list: $SNP_LIST ($(wc -l < "$SNP_LIST") SNPs)"
  echo ""

  if [[ "$INPUT" == *.bed ]]; then
    echo "=== Plink to Plink: $PREFIX ==="
    if [[ $(wc -l < "$SNP_LIST") -gt 20000 ]]; then
      echo "Capping SNP list to 20k random variants (was $(wc -l < "$SNP_LIST"))..."
      shuf -n 20000 "$SNP_LIST" > snp_list_capped.txt
      SNP_LIST="snp_list_capped.txt"
    fi
    INPUT_FLAGS="--bfile ${INPUT%.bed}"
  else
    echo "=== VCF to Plink: $PREFIX ==="
    INPUT_FLAGS="--vcf $INPUT --double-id --max-alleles 2"
  fi

  plink2 \
    $INPUT_FLAGS \
    --extract "$SNP_LIST" \
    --autosome \
    --maj-ref force \
    --make-bed \
    --out "$PREFIX" \
    --threads $THREADS \
    --memory ~{memory_gb * 1024} \
    --allow-extra-chr

  echo ""
  echo "Renaming IIDs to ${PREFIX}_OLDIID..."
  awk -v prefix="$PREFIX" 'BEGIN{OFS="\t"} {print $1, $2, $1, prefix "_" $2}' "${PREFIX}.fam" > id_mapping.txt
  plink2 \
    --bfile "$PREFIX" \
    --update-ids id_mapping.txt \
    --make-just-fam \
    --out "$PREFIX"
  rm -f id_mapping.txt

  echo ""
  echo "Done."
  echo "  SNPs:    $(wc -l < "${PREFIX}.bim")"
  echo "  Samples: $(wc -l < "${PREFIX}.fam")"
  >>>

  output {
    Array[File] plink_data = ["~{prefix}.bed", "~{prefix}.bim", "~{prefix}.fam"]
    File        plink_bim  = "~{prefix}.bim"
  }

  runtime {
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_size} HDD"
    cpu: 16
  }
}
