#!/bin/bash
if [[ $# -ne 1 ]]
then
	echo "  Description: extract LARMORD features"
	echo "  Usage: $0 < PDBIDs>"
	echo "  Usage: $0 2LX1"
else	
	PDB=$1
	
	# Setup Environment
	source ~/.bashrc
	larmorDHome=~/GitSoftware/LarmorD_New/
	
	# Get native ensemble
	mkdir -p ${PDB}	
	split_pdb.bash ${PDB} ${PDB}/model 0 &> ${PDB}/split.out
	models=`grep contains ${PDB}/split.out | awk '{print $3}'`	
	models=`seq 1 $models`	

  # Get data
	data=${dataHome}/measured_shifts_corrected_${PDB}.dat  
	cat $data | sort -n -k2 > ${PDB}/cs_1.txt
	data=${PDB}/cs_1.txt
  
  # Loop over NMR bundle and extract structure features
  for m in $models
  do
  	id=`echo "${PDB}_${m}"`
  	[[ $m == 1 ]] && ${larmorDHome}/bin/larmord_extractor -rna -ring -cutoffRing 10.0 -identification ${id} -csfile ${data} -beta -3 -cutoff 9999 ${PDB}/model_${m}.pdb  	
  	[[ $m != 1 ]] && ${larmorDHome}/bin/larmord_extractor -rna -ring -cutoffRing 10.0 -identification ${id} -csfile ${data} -beta -3 -cutoff 9999 ${PDB}/model_${m}.pdb  | grep -v frame	
  done
fi
