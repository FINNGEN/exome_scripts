version development

workflow exome_convert {
  input{
    File chrom_file_list
    String name
    String docker
    Boolean test
    File mapping
    File variants
    
  }

  Int disk_factor = 5
  Array[Array[String]] chrom_list = read_tsv(chrom_file_list)
  # subset to test scenario
  Array[Array[String]] final_list = chrom_list
  #Array[Array[String]] final_list =  [chrom_list[0]] 
  scatter (elem in final_list){
    call chrom_convert {
      input :
      name = name,
      test=test,
      mapping = mapping,
      docker = docker,
      chrom = elem[0],
      cFile = elem[1],
      variants = variants,
      disk_factor = disk_factor
    }
    call bgen_plink {
      input :
      docker = docker,
      vcf = chrom_convert.vcf,
      chrom = elem[0],
      name = name,
      variants = variants,
      disk_factor = disk_factor +2
    }
  }
}



task bgen_plink {
  input {
    String docker
    File vcf
    String bargs
    String pargs
    File variants
    String chrom
    String name
    Int disk_factor
    File min_pheno
  }
  File tbi = vcf +'.tbi'
  Int vcf_size = ceil(size(vcf,"GB"))
  Int disk_size = disk_factor *vcf_size  + 50
  Int cpu =  8 + 8*(vcf_size/30) + 16*(vcf_size/100)
  String name_chrom =  name + "_" + chrom

  command <<<
  zcat -f  ~{variants} | sed 's/chrX/chr23/g' | sed 's/chrY/chr24/g' | grep chr~{chrom}_ | sed 's/chr23/chrX/g' | sed 's/chr24/chrY/g'   > ./variants.txt
  head ./variants.txt

  # NEEDED FOR CHR23/PAR
  # get SEX COL ID
  COL=$(zcat ~{min_pheno} | head -n1 | tr '\t' '\n' | grep -n SEX | cut -d : -f 1)
  zcat ~{min_pheno} | cut -f 1,$COL | sed -E 1d | sort > sex_info.txt
  bcftools query -l ~{vcf} > samples.txt
  # BUILD SEX UPDATE
  echo -e "FID\tIID\tSEX" > sex_update.txt  && join <(sort samples.txt) sex_info.txt | awk '{print $1"\t"$1"\t"$2}' >> sex_update.txt

  python3 /Scripts/annotate.py  --cFile ~{vcf} --tbiFile ~{tbi}  --oPath "/cromwell_root/"  --vcf-variants ./variants.txt  --name ~{name_chrom}   --split   -pb  --bargs ~{bargs} --pargs ~{pargs} --pconvargs  ' --update-sex sex_update.txt  --split-par b38 '  | tee  chrom_convert_~{name_chrom}.log        
  df -h >> chrom_convert_~{name_chrom}.log
  >>>
  runtime {
    docker: "${docker}"
    cpu: "${cpu}"
    disks: "local-disk ${disk_size} HDD"
    bootDiskSizeGb: 20
    memory: "${cpu}  GB"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: 0
  }

  output {
    # BGEN
    File bgen = "/cromwell_root/${name_chrom}/${name_chrom}.bgen"
    File bgen_sample = "/cromwell_root/${name_chrom}/${name_chrom}.bgen.sample"
    File bgen_index = "/cromwell_root/${name_chrom}/${name_chrom}.bgen.bgi"
    # PLINK
    File bed = "/cromwell_root/${name_chrom}/${name_chrom}.bed"
    File bim = "/cromwell_root/${name_chrom}/${name_chrom}.bim"
    File fam = "/cromwell_root/${name_chrom}/${name_chrom}.fam"
    # LOGS
    File bgen_plink_convert_log  = "chrom_convert_${name_chrom}.log"
    Array[File] logs = glob("/cromwell_root/${name_chrom}/logs/*")

  }
}

 

task chrom_convert {
  input {
    Boolean test
    File mapping
    String chrom
    File cFile
    File variants
    String docker
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
  zcat -f  ~{variants} | sed 's/chrX/chr23/g' | sed 's/chrY/chr24/g' | grep chr~{chrom}_ | sed 's/chr23/chrX/g' | sed 's/chr24/chrY/g'   > ./variants.txt
  head ./variants.txt

  # SAMPLE LIST
  cut -f1  ~{mapping} ~{if test then " | head -n 100" else ""} > ./samples.txt
  wc -l samples.txt
  
  # SUBSAMPLE AND ADD FILTERS
  python3 /Scripts/annotate.py  --cFile ~{cFile} --tbiFile ~{tbiFile}  --oPath "/cromwell_root/"  --vcf-variants ./variants.txt  --name tmp   --split    -v  --samples ./samples.txt   --check-vcf --cleanup  --vargs ~{vargs} | tee  chrom_convert_~{name_chrom}.log        
  df -h >> chrom_convert_~{name_chrom}.log
  
  # REHEADER
  echo "REHEADER"
  bcftools reheader ./tmp/tmp.vcf.gz -s ~{mapping} > ~{name_chrom}.vcf.gz
  tabix ~{name_chrom}.vcf.gz
  df -h >> chrom_convert_~{name_chrom}.log
  >>>
  
  runtime {
    docker: "~{docker}"
    cpu: "~{cpu}"
    disks: "local-disk ~{disk_size} HDD"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    memory:  "~{cpu} GB"
    preemptible: 1
  }
  output {
    File chrom_convert_log  = "chrom_convert_${name_chrom}.log"
    File vcf     = "/cromwell_root/~{name_chrom}.vcf.gz"
    File tbi      = "/cromwell_root/~{name_chrom}.vcf.gz.tbi"
  }
}

