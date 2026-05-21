version 1.0

workflow exome_duplicates {
  input {
    Array[Array[String]] vcf_pairs
    File   plink_bed
    String plink_prefix
    File?  snp_list
  }

  File plink_bim = sub(plink_bed, "\\.bed$", ".bim")
  File plink_fam = sub(plink_bed, "\\.bed$", ".fam")

  Array[File] plink_input_files =  [plink_bed, plink_bim, plink_fam, plink_afreq]

  if (!defined(snp_list)) {
    call ExtractSnpsFromBim {
      input:
        bim = plink_bim
    }
  }

  File actual_snp_list = select_first([snp_list, ExtractSnpsFromBim.snp_list])
  
  call ConvertToPlink as PlinkFilter {
    input:
      prefix      = plink_prefix,
      input_files = plink_input_files,
      snp_list    = actual_snp_list
  }

  scatter (pair in vcf_pairs) {
    call ConvertToPlink as VcfToPlink {
      input:
        prefix      = pair[0],
        input_files = [pair[1], pair[2]],
        snp_list    = actual_snp_list
    }

    call RunKinship {
      input:
        vcf_plink = VcfToPlink.plink_data,
        ref_plink = PlinkFilter.plink_data,
        prefix    = pair[0] + "_king"
    }
  }

  
  output {
    Array[Array[File]] exome_plink   = VcfToPlink.plink_data
    Array[File]        kinship_con   = RunKinship.con_file
  }
}

task ExtractSnpsFromBim {
  input {
    File bim
  }

  Int disk_size = ceil(size(bim, 'GB')) + 5

  command <<<
  set -euo
  cut -f2 ~{bim} > snp_list.txt
  echo "Extracted $(wc -l < snp_list.txt) SNPs from bim file"
  >>>

  output {
    File snp_list = "snp_list.txt"
  }

  runtime {
    memory: "8G"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: 1
  }
}

task RunKinship {
  input {
    Array[File] vcf_plink
    Array[File] ref_plink
    String      prefix
    Int         cpu = 4
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
  king -b "${VCF_BED}","${REF_BED}" --duplicate --prefix "$OUTPUT"

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
    memory: "16 GB"
    disks: "local-disk ~{disk_size} HDD"
    cpu: cpu
    preemptible: 1
  }
}

task ConvertToPlink {
  input {
    String      prefix
    Array[File] input_files
    File        snp_list
  }

  Boolean is_plink  = basename(input_files[0], ".bed") != basename(input_files[0])
  Int     disk_size = ceil(size(input_files[0], 'GB') * 3) + 20

  command <<<
  set -euo
  THREADS=$(nproc)
  PREFIX="~{prefix}"
  SNP_LIST="~{snp_list}"
  INPUT_FILES=(~{sep=" " input_files})
  INPUT="${INPUT_FILES[0]}"

  echo "SNP list: $SNP_LIST ($(wc -l < "$SNP_LIST") SNPs)"
  echo ""

  if [[ "$INPUT" == *.bed ]]; then
    echo "=== Plink to Plink: $PREFIX ==="
    INPUT_FLAGS="--bfile ${INPUT%.bed} --read-freq ${INPUT_FILES[3]}"
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
    Array[File] plink_data = ["~{prefix}.bed", "~{prefix}.bim","~{prefix}.fam"]
  }

  runtime {
    memory: "16 GB"
    disks: "local-disk ~{disk_size} HDD"
    cpu: 16
    preemptible: 1
  }
}
