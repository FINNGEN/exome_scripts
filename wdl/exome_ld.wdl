version development


workflow exome_ld {
  input {
    String docker
    File exome_vcf_files
    File fg_vcf_files
  }

  Array[Array[String]] inputs = read_tsv(exome_vcf_files)
  Map[String,File] vcf_map = read_map(fg_vcf_files)
  scatter (elem in inputs) {
    String chrom = elem[0]
    String exome_plink_root = sub(elem[1],".bed","")
    # subset FG data only to exome samples (after renaming exome data to FG ids)
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
  plink2 --bfile ~{sub(exome_files[0],".bed","")} --pmerge ~{sub(fg_files[0],".bed","")} --make-bed --out ~{out_root}
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
  }
  File exome_fam = exome_plink_root + ".fam"
  Int disk_size = ceil(size(vcf,"GB"))*4
  Int cpus = 8
  String out_root = "finngen_exome_" + chrom
  command <<<
  plink2 --vcf ~{vcf}  ~{pargs}  --keep ~{exome_fam} --make-bed --out ~{out_root}
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
