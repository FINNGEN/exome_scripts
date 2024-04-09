#!/bin/bash

OUT_DIR=$1
KIN=$2
KIN_KIN=$KIN"0"

DUP_LIST=${3-"/mnt/disks/data/samples/exclusions/finngen_R12_duplicate_list.txt"}
EX_SAMPLES=${4-"/mnt/disks/data/exome/plink/kinship/exome_merged.fam"}	    
FG_SAMPLES=${5-"/home/pete/nfs/fgfactory_pass_samples/R12_fgfactory_pass_samples.txt"}

mkdir -p $OUT_DIR
FILTERED_KIN="$OUT_DIR/merged.kin"
echo $FILTERED_KIN

# WRITE ONLY COUPLES ACROSS TWO DATA SETS
cut -f 2,3,15 $KIN | grep "QRY_" | grep "REF_"  > $FILTERED_KIN.tmp &&   cut -f 2,4,14 $KIN_KIN | grep "QRY_" | grep "REF_" >> $FILTERED_KIN.tmp
# fix order so we have on one side FG and the other EX
paste  \
    <( cat $FILTERED_KIN.tmp |   grep -oh "REF_FG\w*"  | cut  -d "_" -f2) \
    <( cat $FILTERED_KIN.tmp |   grep -oh "REF_FG\w*" ) \
    <( cat $FILTERED_KIN.tmp |   grep -oh "QRY_FG\w*" | cut -d "_" -f2) \
    <( cat $FILTERED_KIN.tmp |   grep -oh "QRY_FG\w*") \
    <( cat $FILTERED_KIN.tmp |   cut -f 3 ) \
					  > $FILTERED_KIN
 
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
 
