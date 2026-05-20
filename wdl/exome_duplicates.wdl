version 1.0

workflow exome_duplicates {
  input {
    Array[Pair[String, String]] vcf_pairs
    File   plink_bed
    String plink_prefix
    File?  snp_list
    File?  plink_afreq
  }

  File plink_bim = sub(plink_bed, "\\.bed$", ".bim")
  File plink_fam = sub(plink_bed, "\\.bed$", ".fam")

  if (!defined(snp_list)) {
    call ExtractSnpsFromBim {
      input:
        bim = plink_bim
    }
  }

  File actual_snp_list = select_first([snp_list, ExtractSnpsFromBim.snp_list])

  scatter (pair in vcf_pairs) {
    call ConvertToPlink as VcfToPlink {
      input:
        prefix     = pair.left,
        input_file = pair.right,
        snp_list   = actual_snp_list
    }
  }

  call ConvertToPlink as PlinkFilter {
    input:
      prefix        = plink_prefix,
      input_file    = plink_bed,
      snp_list      = actual_snp_list,
      sidecar_bim   = plink_bim,
      sidecar_fam   = plink_fam,
      sidecar_afreq = plink_afreq
  }

  output {
    Array[Array[File]] exome_plink   = VcfToPlink.plink_data
    Array[File] fg_plink             = PlinkFilter.plink_data
    File             snp_list_used   = actual_snp_list
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

task ConvertToPlink {
  input {
    String prefix
    File   input_file
    File   snp_list
    File?  sidecar_bim
    File?  sidecar_fam
    File?  sidecar_afreq
  }

  Boolean is_plink  = basename(input_file, ".bed") != basename(input_file)
  Int     disk_size = ceil(size(input_file, 'GB') * 3) + 20

  command <<<
  set -euo
  THREADS=$(nproc)
  PREFIX="~{prefix}"
  INPUT="~{input_file}"
  SNP_LIST="~{snp_list}"
  AFREQ="~{if defined(sidecar_afreq) then select_first([sidecar_afreq]) else ""}"

  echo "SNP list: $SNP_LIST ($(wc -l < "$SNP_LIST") SNPs)"
  echo ""

  if [[ "$INPUT" == *.bed ]]; then
    echo "=== Plink to Plink: $PREFIX ==="
    INPUT_FLAGS="--bfile ${INPUT%.bed}${AFREQ:+ --read-freq $AFREQ}"
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
