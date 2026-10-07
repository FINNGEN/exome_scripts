version development


workflow assign_ids {

  input {
    String docker
    File kin
    File dup_list
    File ex_samples
    File fg_samples
  }

  call kinship {
    input :
    docker= docker,
    kin= kin,
    dup_list= dup_list,
    ex_samples= ex_samples,
    fg_samples= fg_samples
  }
  call assign {
    input :
    docker= docker,
    dup_list = dup_list,
    exome_fam = ex_samples,
    matching = kinship.matching,
    tentative = kinship.tentative,
    merged_kin = kinship.merged_kin
  }

}

task assign {
  input {
    String docker
    File exome_fam
    File merged_kin
    File tentative
    File matching
    File dup_list
  }
  
  command <<<

  mkdir ./results/
  cat << "__EOF__" > assign.py
  import sys,os

  out_dir=sys.argv[1]
  # INPUTS
  exome_fam=sys.argv[2]
  kinship_file=sys.argv[3]
  tentative_matching_list=sys.argv[4]
  matching_file=sys.argv[5]
  dup_list=sys.argv[6]

  #OUTPUTS
  rej_file=os.path.join(out_dir,"problematic_relations.txt")
  out_map=os.path.join(out_dir,"exome_fg_id_mapping.txt" )
  out_rej=os.path.join(out_dir,"exome_fg_id_rejected_issues.txt" )

  tent_map = {}
  with open(tentative_matching_list) as i:
      for line in i:
          ids = line.strip().split()
          id1,id2 = [ids[0],ids[2]]
          tent_map[id1] = id2
          tent_map[id2] = id1



  with open(matching_file) as i: matching_samples=set([line.strip() for line in i])

  # GOTTA CHANGE LOGIC
  # FIRST I PASS AND ASSIGN ALL THE MATCHING DUPLICATES AND TENTATIVE DUPLICATES. THEN I DO ANOTHER PASS AND TRY TO SORT OUT THE REDUNDANT DUPLICATES
  #now i dothe real assigning
  id_mapping = {}
  with open(matching_file) as i:
      for elem in set([line.strip() for line in i]):
          id_mapping[elem] = elem
  for elem in tent_map:
      id_mapping[elem] = tent_map[elem]

  print(f"{len(id_mapping)} have ok mapping")


  # let's filter the known duplicate list only to samples that are needed in the analysis
  # here i return the list of FG ids from the kinship data who have a Dup/MZ (i.e. will need to be investigated)
  with open(kinship_file) as i:
      fg_ex_dups = []
      for line in i:
          fg,_,ex_root,ex_full,kinship = line.strip().split()
          if kinship == "Dup/MZ":
              fg_ex_dups.append(fg)

  fg_ex_dups = set(fg_ex_dups)
  print(len(fg_ex_dups))


  # known duplicate mapping
  dup_map = {}
  with open(dup_list) as i:
      next(i)
      for line in i:
          ids = line.strip().split()
          # check if any of the ids is in FG
          check = [(elem in fg_ex_dups) for elem in ids]
          if any(check):
              for id in ids:
                  matches = [elem for elem in ids if elem != id]
                  dup_map[id] = matches

  issue_dict = {}
  # here i do a pass where i label ids that should be DUPs (same id) but don't match genetically
  with open(kinship_file) as i,open(out_rej,'wt') as rej:
      for line in i:
          fg,_,ex_root,ex_full,kinship = line.strip().split()
          ex_full = ex_full.split("QRY_")[1]
          # so where we wanna skip who we fixed already
          if kinship != "Dup/MZ" and fg == ex_root and fg not in tent_map:
              issue_dict[ex_full] = "SAME_ID_BUT_MISMATCH"


  # so now we have already assigned the "ovbious ones" and we need to find if there's a way to map duplicates via known duplicate lists
  with open(kinship_file) as i :
      for line in i:
          fg,_,ex_root,ex_full,kinship = line.strip().split()
          ex_full = ex_full.split("QRY_")[1]
          # so where we wanna skip who we fixed already
          if kinship == "Dup/MZ" and ex_full not in id_mapping and ex_full not in issue_dict:
              assert fg!=ex_root
              # here we have a duplicat where
              # 1) there is no matching IDs (took care of it before)
              # gotta find, if it exists the intersection
              candidate_map = ""
              if fg in dup_map:
                  if ex_root in dup_map[fg]:
                      candidate_map = fg
              # so now i check if it's the other way
              elif ex_root in dup_map:
                  if fg in dup_map[ex_root]:
                      candidate_map = fg
              else:
                  issue_dict[ex_full] = "DUP_BUT_LACKING_ID_MAPPING"
              if candidate_map:
                  # now i need to check if the element has already been added and in case add a v2/v3 label to keep track!
                  counts = len([elem for elem in id_mapping if id_mapping[elem].split("_")[0] == candidate_map])
                  label = f"_exdup{counts+1}" if counts >0  else ""
                  # now i update the mapping, adding label if needed        
                  id_mapping[ex_full] = candidate_map + label
              else:
                  issue_dict[ex_full] = "DUP_BUT_LACKING_ID_MAPPING"
  
  # write mapping to file (need to include dup)
  with open(exome_fam) as i: exome_samples = [elem.strip().split()[0] for elem in i.readlines()]
  with open(out_map,'wt') as acc,open(out_rej,'wt') as rej:
      for elem in exome_samples:
          # log problematic samples
          if elem in issue_dict:
              rej.write("\t".join([elem,issue_dict[elem]]) + "\n")
          elif elem in id_mapping:
              acc.write("\t".join([elem,elem,id_mapping[elem].split("_")[0],id_mapping[elem]]) + '\n')
          else:
              rej.write("\t".join([elem,"NONE_missing"]) + "\n")

  def lines(f):
      with open(f, 'r') as fp: return len(fp.readlines())

  print(lines(exome_fam))
  print(lines(out_map) + lines(out_rej))
  
  __EOF__
  
  python3 assign.py ./results/ ~{exome_fam} ~{merged_kin} ~{tentative} ~{matching} ~{dup_list}
  ls ./results/
  sort -u -k 3,3 ./results/exome_fg_id_mapping.txt | cut -f1,3> final_mapping.txt
  >>>
  output {
    File mapping = "./final_mapping.txt"
    Array[File] misc = glob("./results/*")
  }
  runtime {
    docker : "~{docker}"
    disks: "local-disk 1  HDD"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: 1
  }
}

task kinship {
  
  input {
    String docker
    File kin
    File dup_list
    File ex_samples
    File fg_samples
  }
  File kin0 = kin + "0"

  command <<<
  OUT_DIR="./kinship/"
  KIN=~{kin}
  KIN_KIN=$KIN"0"

  DUP_LIST=~{dup_list}
  EX_SAMPLES=~{ex_samples}
  FG_SAMPLES=~{fg_samples}

  mkdir -p $OUT_DIR
  FILTERED_KIN="$OUT_DIR/merged.kin"
  echo $FILTERED_KIN

  # WRITE ONLY COUPLES ACROSS TWO DATA SETS
  cut -f 2,3,15 $KIN | grep "QRY_" | grep "REF_"  > $FILTERED_KIN.tmp &&   cut -f 2,4,14 $KIN_KIN | grep "QRY_" | grep "REF_" >> $FILTERED_KIN.tmp
  # fix order so we have on one side FG and the other EX
  paste     \
      <( cat $FILTERED_KIN.tmp |   grep -oh "REF_FG\w*"  | cut  -d "_" -f2)  <( cat $FILTERED_KIN.tmp |   grep -oh "REF_FG\w*" )      <( cat $FILTERED_KIN.tmp |   grep -oh "QRY_FG\w*" | cut -d "_" -f2)        <( cat $FILTERED_KIN.tmp |   grep -oh "QRY_FG\w*")  <( cat $FILTERED_KIN.tmp |   cut -f 3 )     > $FILTERED_KIN

  echo "number of cross data set relations"
  cat $FILTERED_KIN | wc -l

  # first let's address what is alright,
  #that is the guys with same ID and DUP
  OK_DUP=$OUT_DIR/"matching.txt"
  echo "matching ids"
  awk 'BEGIN{OFS="\t"} $1==$3' $FILTERED_KIN | grep "Dup/MZ" | cut -f 4 |sed 's/QRY_//g' | sort > $OK_DUP
  cat $OK_DUP | wc -l
  echo "Sanity check with input fam"
  comm -12 <(cut -f2 $EX_SAMPLES | sort | uniq ) $OK_DUP | wc -l

  #DUPLICATES WITH WRONG SAMPLES

  # LET'S CHECK FOR DUPLICATES
  FG_DUP=$OUT_DIR/"fg_duplicates.txt"
  cat $DUP_LIST | tr "\t" "\n" | sort | uniq > $FG_DUP
  echo "checking Dup/MZ with non matching ID but also not in OK list as there might be twins as inputs"
  NOT_OK_DUP=$OUT_DIR/"wrong_matching.txt"
  awk 'BEGIN{OFS="\t"} $1!=$3' $FILTERED_KIN | grep -vwf $OK_DUP | grep "Dup/MZ"  >$NOT_OK_DUP
  cat $NOT_OK_DUP | wc -l

  # these are samples that have 2 twins in FG but neither of them matches. It turns out that these are scenarios where they are in the FG duplicate list of a known twin. That is, the exome id is the duplicate ID that does not appear in FG of an actual twin (thus two matches, one being the duplicate, the other the twin). There are also cases when the twin couplee has *BOTH* wrong sides of the duplicate list
  PROB_TWINS=$OUT_DIR/"problematic_twins.txt"
  cut -f4 $NOT_OK_DUP | sort | uniq -D | sort | uniq -c | sed 's/QRY_//g'   | awk '{print $2}' > $PROB_TWINS


  # GET LIST OF SAMPLES THAT ARE NOT IN DUPLICATES (EITHER FG OR EXOME)
  cut -f1,3 $NOT_OK_DUP | tr "\t" "\n" | grep -vwf $FG_DUP | sort | uniq  > $OUT_DIR/non_dup.txt
  echo "samples that do not have a known duplicate"
  cat $OUT_DIR/non_dup.txt | wc -l
  #cat $NOT_OK_DUP | grep -vwf $OUT_DIR/non_dup.txt | cut -f4  | sort | uniq |wc -l
  # WE NEED COUPLE WHERE BOTH SAMPLES ARE NOT IN DUP LIST

  # samples that cannot be explained between duplicates and problematic twins
  UNEXPLAINED_WRONG=$OUT_DIR/"wrong_matching_unexplained.txt"
  cat $NOT_OK_DUP | grep -wvf $PROB_TWINS | grep -wf $OUT_DIR/non_dup.txt | cut -f4 | sed 's/QRY_//g' > $UNEXPLAINED_WRONG
  echo "samples that have a twin and there is no explanation through duplicate list"
  cat $UNEXPLAINED_WRONG | wc -l

  #PROBABLE MISMATCHES
  # these are cases where the misssed couple appears also in the other direction
  TENT_MISMATCH=$OUT_DIR/"tentative_mismatches.txt"
  cat $NOT_OK_DUP | grep -wf <(cat $NOT_OK_DUP | grep -wf $UNEXPLAINED_WRONG  | cut -f1,3 | tr "\t" "\n" | sort | uniq -c | awk '{print $1"\t"$2}' | awk '$1>1' | cut -f2 )> $TENT_MISMATCH

  MIS=$OUT_DIR/"missing.txt"
  comm -23 <(cut -f2 $EX_SAMPLES | sort ) $OK_DUP | comm -23 -  <( cut -f4 $NOT_OK_DUP | sed 's/QRY_//g' | sort)  > $MIS

  N_MATCH=$(cat $OK_DUP | wc -l)
  N_NO_MATCH=$(cat $NOT_OK_DUP| cut -f4 $NOT_OK_DUP |sed 's/QRY_//g' | sort | uniq | wc -l)
  N_MISS=$(cat $MIS | wc -l)
  num=$((N_MATCH + N_NO_MATCH + N_MISS))
  echo $num
  cat $EX_SAMPLES | wc -l

  #OVERLAPPING_MISSIN
  MIS_OV=$OUT_DIR/"missing_overlap.txt"
  MIS_NO_OV=$OUT_DIR/"missing_no_overlap.txt"

  comm -12 <(sort $MIS) <(cut -d ":" -f6 $FG_SAMPLES |sort) > $MIS_OV
  comm -23 <(sort $MIS) <(cut -d ":" -f6 $FG_SAMPLES |sort) > $MIS_NO_OV

  rm $FILTERED_KIN.tmp  $FG_DUP $OUT_DIR/non_dup.txt

  ls ./kinship/*
  >>>
  runtime {
    docker : "~{docker}"
    disks: "local-disk 1  HDD"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: 1
    
  }

  output {
    File merged_kin = "./kinship/merged.kin"
    File matching = "./kinship/matching.txt"
    File tentative = "./kinship/tentative_mismatches.txt"
    Array[File] all = glob("./kinship/*")
  }

}
 
