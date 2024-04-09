
import sys,os

data_dir=sys.argv[1]
# INPUTS
exome_fam=sys.argv[2]
kinship_file=os.path.join(data_dir,"merged.kin")
tentative_matching_list= os.path.join(data_dir,"tentative_mismatches.txt")
matching_file=os.path.join(data_dir,"matching.txt")
dup_list="/mnt/disks/data/samples/exclusions/finngen_R12_duplicate_list.txt"

#OUTPUTS
rej_file=os.path.join(data_dir,"problematic_relations.txt")
out_map=os.path.join(data_dir,"exome_fg_id_mapping.txt" )
out_rej=os.path.join(data_dir,"exome_fg_id_rejected_issues.txt" )


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
        if kinship != "Dup/MZ" and fg == ex_root:
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
