version development


workflow exome_ld {
  input {
    String docker
    File exome_beds
    File finngen_vcfs
  }
  # READ IN ALL EXOME BEDS
  Array[Array[String]] inputs = read_tsv(exome_beds)
  # MAP TO ALL VCFS, SO I CAN WORK WITH SUBSET IF NEEDED
  Map[String,File] vcf_map = read_map(finngen_vcfs)

  Array[Array[String]] final_inputs = [inputs[22]]
  scatter (elem in inputs) {
    String chrom = elem[0]
    String exome_plink_root = sub(elem[1],".bed","")
    # subset FG data only to exome samples (exome data must have been renamed)
    call filter_fg_plink{
      input : docker = docker,chrom = chrom,vcf = vcf_map[chrom],exome_plink_root = exome_plink_root
    }
    # filter out FG variants from exome
    call filter_exome_variants {
      input : docker = docker,chrom=chrom,exome_plink_root = exome_plink_root,fg_bim = filter_fg_plink.bim
    }
    Array[File] exome_files = [filter_exome_variants.bed,filter_exome_variants.bim,filter_exome_variants.fam]
    Array[File] fg_files    = [ filter_fg_plink.bed,filter_fg_plink.bim,filter_fg_plink.fam]
    # merge into single dataset
    call merge_files {
      input:docker = docker,exome_files=exome_files,fg_files=fg_files,chrom=chrom
    }
    call ld {
      input: docker = docker,chrom=chrom,plink_files = [merge_files.bed,merge_files.bim,merge_files.fam],exome_bim = exome_files[1],finngen_bim = fg_files[1]
    }
  }
}


task ld {
  input {
    String docker
    Array[File] plink_files
    File exome_bim
    File finngen_bim
    String ld_params
    String chrom
  }

  Int disk_size = ceil(size(plink_files,"GB"))*2
  Int cpus = 8
  String out_root = "exome_finngen_ld_" + chrom
  command <<<
  # calculate R2 only for exome variants
  plink2 --bfile ~{sub(plink_files[0],'.bed','')} --r2-unphased ~{ld_params} --ld-snp-list ~{exome_bim} --out ld
  echo -e "EXOME_SNP\tFINNGEN_SNP\tR2" > ~{out_root}.ld
  # join over FG variants
  join -1 2 <(cut -f 3,6,7 ld.vcor | sort -k2 ) <(cut -f2 ~{finngen_bim} | sort ) | awk '{ print $2"\t"$1"\t"$3}' | sort -k1 >> ~{out_root}.ld
  
  >>>

  runtime {
    docker: "${docker}"
    cpu: "${cpus}"
    disks: "local-disk " + "${disk_size}" + " HDD"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    memory:  "${cpus} GB"
    preemptible: 1
  }
  output {
    File ld         = "~{out_root}.ld"
  }
}

task merge_files {
  input {
    String docker
    Array[File] exome_files
    Array[File] fg_files
    String chrom
  }
  
  String out_root = "exome_fg_merged_" + chrom
  Int disk_size = ceil(size(exome_files,"GB") + size(fg_files,"GB"))*4+10
  Int cpus = 8
  command <<<
  plink --bfile ~{sub(exome_files[0],".bed","")} --bmerge ~{sub(fg_files[0],".bed","")} --make-bed --out ~{out_root} --allow-extra-chr
  plink2 --bfile ~{out_root} --make-just-fam --out tmp
  mv tmp.fam ~{out_root}.fam
  >>>
  runtime {
    docker: "${docker}"
    cpu: "${cpus}"
    disks: "local-disk " + "${disk_size}" + " HDD"
    bootDiskSizeGb: 20
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    memory:  "${cpus} GB"
    preemptible: 1
  }
  output {
    File bed         = "~{out_root}.bed"
    File bim         = "~{out_root}.bim"
    File fam         = "~{out_root}.fam"     
  }
}

task filter_exome_variants {

  input {
    String docker
    String chrom
    String exome_plink_root
    File fg_bim
  }

  Array[File] exome_plink_files = [exome_plink_root + ".bed",exome_plink_root +".bim", exome_plink_root + ".fam"]
  Int disk_size = ceil(size(exome_plink_files[0],"GB"))*4 + 10
  Int cpus = 8
  String out_root = "exome_no_fg_variants_" + chrom
  command <<<
  plink2 --bfile ~{sub(exome_plink_files[0],".bed","")} --exclude ~{fg_bim} --make-bed --out ~{out_root}
  >>>

  runtime {
    docker: "${docker}"
    cpu: "${cpus}"
    disks: "local-disk " + "${disk_size}" + " HDD"
    bootDiskSizeGb: 20
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    memory:  "${cpus} GB"
    preemptible: 1
  }
  output {
    File bed         = "~{out_root}.bed"
    File bim         = "~{out_root}.bim"
    File fam         = "~{out_root}.fam"     
  }
}

task filter_fg_plink {
  input {
    String docker
    String chrom
    File vcf
    String exome_plink_root
    String pargs
    File min_pheno
  }
  File exome_fam = exome_plink_root + ".fam"
  Int disk_size = ceil(size(vcf,"GB"))*4
  Int cpus = 8
  String out_root = "finngen_exome_samples_" + chrom
  command <<<
  COL=$(zcat ~{min_pheno} | head -n1 | tr '\t' '\n' | grep -n SEX | cut -d : -f 1)
  zcat ~{min_pheno} | cut -f 1,$COL | sed -E 1d | sort > sex_info.txt
  bcftools query -l ~{vcf} > samples.txt
  # BUILD SEX UPDATE
  echo -e "FID\tIID\tSEX" > sex_update.txt  && join <(sort samples.txt) sex_info.txt | awk '{print $1"\t"$1"\t"$2}' >> sex_update.txt
  
  plink2 --vcf ~{vcf}  ~{pargs}  --keep ~{exome_fam} --make-bed --out ~{out_root} --update-sex sex_update.txt  --split-par b38
  >>>

  runtime {
    docker: "${docker}"
    cpu: "${cpus}"
    disks: "local-disk " + "${disk_size}" + " HDD"
    bootDiskSizeGb: 20
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    memory:  "${cpus} GB"
    preemptible: 1
  }
  output {
    File bed         = "~{out_root}.bed"
    File bim         = "~{out_root}.bim"
    File fam         = "~{out_root}.fam"     
  }

}
