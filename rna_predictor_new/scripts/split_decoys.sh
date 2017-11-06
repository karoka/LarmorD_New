#!/bin/bash

# script to split decoy files
source ~/.bashrc	
rnas="1L1W 1LC6 1LDZ 1NC0 1OW9 1PJY 1R7W 1R7Z 1UUU 1YSV 2FDT 2GM0 2JXQ 2JXS 2K3Z 2K41 2KOC 2KYD 2LBJ 2LBL 2LDL 2LDT 2LHP 2LI4 2LK3 2LP9 2LPA 2LU0 2LUB 2LV0 2RN1 2Y95 4A4S 4A4T 4A4U 2JWV 2K63 2K64 2K65 2K66 2LPS 2LQZ 2LUN 2LX1 2M12 2M21 2M22 2M23 2M24 2M57 2M8K 2MEQ 2MHI 2MI0 2MIS 1SCL 1XHP 1ZC5 2QH2 2QH3 2QH4 1Z2J 2JWV 2K63 2K64 2K65 2K66 2LPS 2LQZ 2LUN 2LX1 2M12 2M21 2M22 2M23 2M24 2M57 2M8K 2MEQ 2MHI 2MI0 2MIS 1SCL 1XHP 1ZC5 2QH2 2QH3 2QH4"

# loop over rnas and get data
for rna in ${rnas}
do	
	inputname=${rna}/${rna}_decoys_1_nlb_decoys
  mkdir -p decoy_files/${rna}
  
	# Count number of models using catdcd
	models=`catdcd -pdb ${inputname}.pdb | tail -2 | head -1 | awk '{print $3}'`	
	echo "${rna} contains ${models} models"	
	
	# Create sequence of model numbers
	models=`seq 1 ${models}`	
	# Loop over models and write out individual PDBs
	for i in $models
	do
		model=$((i+offset))
		catdcd -o decoy_files/${rna}/decoy_${model}.pdb -otype pdb -s ${inputname}.pdb -stype pdb -first ${i} -last ${i} -pdb ${inputname}.pdb #> /dev/null
	done
done	
