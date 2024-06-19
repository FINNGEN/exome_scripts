version development

workflow exome_annotate {
  input{
    File chrom_file_list
    String name
    Boolean test
  }
  Int cpus =64
  Int mem = 32
  Int disk_factor = 4

  Array[Array[String]] chrom_list = read_tsv(chrom_file_list)
  Array[Array[String]] final_list = chrom_list
  # subset to test scenario
  #Array[Array[String]] final_list =  [chrom_list[22],chrom_list[23]] 
  scatter (elem in final_list){
    call chrom_convert {
      input :
      name = name,
      test=test,
      chrom = elem[0],
      cFile = elem[1],
      cpus = cpus,mem = mem,disk_factor = disk_factor
    }
  }
}

task chrom_convert {
  input {
    Boolean test
    String chrom
    File cFile
    File variants
    String missingness
    String vargs
    String name
    Int disk_factor
    Int mem
    Int cpus
  }

  Int disk_size = disk_factor * ceil(size(cFile,"GB")) + 50

  String name_chrom =  name + "_" + chrom
  File tbiFile = cFile + '.tbi'
  command <<<
   zcat -f  ~{variants} | sed 's/chrX/chr23/g' | sed 's/chrY/chr24/g' | grep chr~{chrom}_ | sed 's/chr23/chrX/g' | sed 's/chr24/chrY/g'   > ./variants.txt
   head ./variants.txt

   bcftools query -l  ~{cFile} ~{if test then " | head -n 10" else ""} > ./samples.txt
   wc -l samples.txt
   
   python3 /Scripts/annotate.py  --cFile ~{cFile} --tbiFile ~{tbiFile}  --oPath "/cromwell_root/"  --vcf-variants ./variants.txt  --name ~{name_chrom}   --split    -v   ~{if test then "--samples ./samples.txt" else ""}   --annotate "" --check-vcf  --set-missingness ~{missingness} --vargs ~{vargs} | tee  chrom_convert_~{name_chrom}.log        
  df -h >> chrom_convert_~{name_chrom}.log
  >>>
  
  runtime {
    cpu: "~{cpus}"
    disks: "local-disk ~{disk_size} HDD"
    memory:  "~{mem} GB"
    preemptible: 1
  }
  output {
    File chrom_convert_log  = "chrom_convert_${name_chrom}.log"
    File vcf     = "/cromwell_root/~{name_chrom}/~{name_chrom}.vcf.gz"
    File tbi      = "/cromwell_root/~{name_chrom}/~{name_chrom}.vcf.gz.tbi"
  }
}

