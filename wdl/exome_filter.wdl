version development

workflow exome_filter {
  input{
    File chrom_file_list
    String name
    String docker
    Boolean test
    File mapping
  }

  Array[Array[String]] chrom_list = read_tsv(chrom_file_list)
  Array[Array[String]] final_list = chrom_list
  # subset to test scenario

  Int cpus = 8
  Int mem = 16
  Int disk_factor = 4
  
  Array[Array[String]] final_list =  [chrom_list[23]] 
  scatter (elem in final_list){
    call chrom_convert {
      input :
      name = name,
      test=test,
      mapping = mapping,
      docker = docker,
      chrom = elem[0],
      cFile = elem[1],
      cpus = cpus,mem=mem,disk_factor = disk_factor
    }
    call bgen_plink {
      input :
      docker = docker,
      vcf = chrom_convert.vcf,
      chrom = elem[0],
      name = name,
      cpus = cpus,mem=mem,disk_factor = disk_factor +1
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
    Int cpu
    Int mem
    Int disk_factor
  }
  File tbi = vcf +'.tbi'
  
  Int disk_size = disk_factor * ceil(size(vcf,"GB")) + 50
  String name_chrom =  name + "_" + chrom

  command <<<
  zcat -f  ~{variants} | sed 's/chrX/chr23/g' | sed 's/chrY/chr24/g' | grep chr~{chrom}_ | sed 's/chr23/chrX/g' | sed 's/chr24/chrY/g'   > ./variants.txt
   head ./variants.txt

   python3 /Scripts/annotate.py  --cFile ~{vcf} --tbiFile ~{tbiFile}  --oPath "/cromwell_root/"  --vcf-variants ./variants.txt  --name ~{name_chrom}   --split   -pb  --bargs ~{bargs} --pargs ~{pargs} | tee  chrom_convert_~{name_chrom}.log        
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
    File out_bgen = "/cromwell_root/${name_chrom}/${name_chrom}.bgen"
    File out_bgen_sample = "/cromwell_root/${name_chrom}/${name_chrom}.bgen.sample"
    File out_bgen_index = "/cromwell_root/${name_chrom}/${name_chrom}.bgen.bgi"
    File out_chrom_convert_log  = "chrom_convert_${name_chrom}.log"
    Array[File] logs = glob("/cromwell_root/${name_chrom}/logs/*")
    File bed = "/cromwell_root/${name_chrom}/${name_chrom}.bed"
    File bim = "/cromwell_root/${name_chrom}/${name_chrom}.bim"
    File fam = "/cromwell_root/${name_chrom}/${name_chrom}.fam"
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
    Int cpu
    Int mem
    Int disk_factor
  }

  Int disk_size = disk_factor * ceil(size(cFile,"GB")) + 50
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
    cpu: "~{cpus}"
    disks: "local-disk ~{disk_size} HDD"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    memory:  "~{mem} GB"
    preemptible: 1
  }
  output {
    File chrom_convert_log  = "chrom_convert_${name_chrom}.log"
    File vcf     = "/cromwell_root/~{name_chrom}.vcf.gz"
    File tbi      = "/cromwell_root/~{name_chrom}.vcf.gz.tbi"
  }
}

