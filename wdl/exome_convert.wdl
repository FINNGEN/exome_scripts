version development

workflow exome_convert {
  input{
    File chrom_file_list
    String name
    Boolean test
    File mapping
  }

  Int disk_factor = 3
  Array[Array[String]] chrom_list = read_tsv(chrom_file_list)
  # subset to test scenario
  Array[Array[String]] final_list = chrom_list
  #Array[Array[String]] final_list =  [chrom_list[0]] 
  scatter (elem in final_list){
    String chrom = elem[0]
    
    call chrom_convert {
      input :
      name = name,
      test=test,
      mapping = mapping,
      chrom = chrom,
      cFile = elem[1],
      disk_factor = disk_factor
    }
    call plink {
      input :
      vcf = chrom_convert.vcf,
      chrom = elem[0],
      name = name,
      disk_factor = disk_factor + 1
    }
    call bgen {
      input :
      vcf = chrom_convert.vcf,
      chrom = elem[0],
      name = name,
      disk_factor = disk_factor
    }
    
  }
}


task plink {
  input {
    File vcf
    String pargs
    String chrom
    String name
    Int disk_factor
    File min_pheno
  }
  File tbi = vcf +'.tbi'
  Int vcf_size = ceil(size(vcf,"GB"))
  Int disk_size = disk_factor *vcf_size  + 50
  Int cpu =  8
  String name_chrom =  name + "_" + chrom
  String par = if chrom=="23" then "--split-par b38" else ""  
  command <<<
  # NEEDED FOR CHR23/PAR
  # get SEX COL ID
  COL=$(zcat ~{min_pheno} | head -n1 | tr '\t' '\n' | grep -n SEX | cut -d : -f 1)
  zcat ~{min_pheno} | cut -f 1,$COL | sed -E 1d | sort > sex_info.txt
  bcftools query -l ~{vcf} > samples.txt
  # BUILD SEX UPDATE
  echo -e "FID\tIID\tSEX" > sex_update.txt  && join <(sort samples.txt) sex_info.txt | awk '{print $1"\t"$1"\t"$2}' >> sex_update.txt
  plink2 --vcf ~{vcf} ~{pargs} --freq --make-bed --out ~{name_chrom} --update-sex sex_update.txt   ~{par}

  >>>
  runtime {
    cpu: "${cpu}"
    disks: "local-disk ${disk_size} HDD"
    memory: "${cpu}  GB"
    preemptible: 1
  }

  output {
    # PLINK
    File bed  = "${name_chrom}.bed"
    File bim  = "${name_chrom}.bim"
    File fam  = "${name_chrom}.fam"
    File freq = "${name_chrom}.afreq"
  }
}

task bgen {
  input {
    File vcf
    String bargs
    String chrom
    String name
    Int disk_factor
  }
  File tbi = vcf +'.tbi'
  Int vcf_size = ceil(size(vcf,"GB"))
  Int disk_size = disk_factor *vcf_size  + 50
  Int cpu =  8
  String name_chrom =  name + "_" + chrom

  command <<<
  qctool -g ~{vcf} ~{bargs} -og ~{name_chrom}.bgen -os ~{name_chrom}.bgen.sample
  bgenix -g  ~{name_chrom}.bgen -clobber -index
  >>>
  runtime {
    cpu: "${cpu}"
    disks: "local-disk ${disk_size} HDD"
    memory: "${cpu}  GB"
    preemptible: 1
  }

  output {
    # BGEN
    File bgen        = "${name_chrom}.bgen"
    File bgen_sample = "${name_chrom}.bgen.bgi"
    File bgen_index  = "${name_chrom}.bgen.sample"
  }
}

 

task chrom_convert {
  input {
    Boolean test
    File mapping
    String chrom
    File cFile
    String vargs
    String name
    Int cpu =  8 + 8*(vcf_size/30) + 16*(vcf_size/100)
    Int disk_factor
  }

  Int vcf_size = ceil(size(cFile,"GB"))
  Int disk_size = disk_factor *vcf_size  + 50
  Int cpu =  8 + 8*(vcf_size/20) + 16*(vcf_size/60) #8/16/32 cpus based on size
  String name_chrom =  name + "_" + chrom
  File tbiFile = cFile + '.tbi'
  command <<<
  # BUILD VARIANT LIST
  # SAMPLE LIST
  cut -f1  ~{mapping} ~{if test then " | head -n 100" else ""} > ./samples.txt
  wc -l samples.txt
  # SUBSAMPLE AND ADD FILTERS
  bcftools view ~{cFile} -S samples.txt -Ou | bcftools view ~{vargs} -Oz -o tmp.vcf.gz
  rm ~{cFile}
  # REHEADER
  echo "REHEADER"
  bcftools reheader tmp.vcf.gz -s ~{mapping} > ~{name_chrom}.vcf.gz
  tabix ~{name_chrom}.vcf.gz
  bcftools index -s ~{name_chrom}.vcf.gz
  >>>
  
  runtime {
    cpu: "~{cpu}"
    disks: "local-disk ~{disk_size} HDD"
    memory:  "~{cpu} GB"
    preemptible: 1
  }
  output {
    File vcf = "/cromwell_root/~{name_chrom}.vcf.gz"
    File tbi = "/cromwell_root/~{name_chrom}.vcf.gz.tbi"
  }
}

