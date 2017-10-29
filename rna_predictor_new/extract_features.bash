#!/bin/bash
if [[ $# -ne 2 ]]
then
	echo "  Description: extract LARMORD features"
	echo "  Usage: $0 < PDBIDs>"
	echo "  Usage: $0 2LX1"
else	
	PDB=$1
	cut_off=$2
	
	# Setup Environment
	source ~/.bashrc
	larmorDHome=~/Documents/LarmorD_New/
	data=${larmorDHome}/rna_predictor/data/measured_shifts_corrected_${PDB}.dat
	
	# go to working directory
	cd ${larmorDHome}/rna_predictor/
	
	# Get number of models
	models=`ls ${PDB}/model_*pdb | grep -v min | wc -l | awk '{print $1}'`
	models=`seq 1 ${models}`	
  
	# Loop over NMR bundle and extract structure features
	for m in $models
	do
		#id=`echo "${PDB}_${m}"`
  		#[[ $m == 1 ]] && ${larmorDHome}/bin/larmord_extractor -rna -ring -cutoffRing 10.0 -identification ${id} -csfile ${data} -beta -3 -cutoff $2 ${PDB}/model_${m}.pdb
  		#[[ $m != 1 ]] && ${larmorDHome}/bin/larmord_extractor -rna -ring -cutoffRing 10.0 -identification ${id} -csfile ${data} -beta -3 -cutoff $2 ${PDB}/model_${m}.pdb  | grep -v frame	
         
        id=`echo "${PDB}_${m}_minimize"`
        [[ $m == "1" ]] && ${larmorDHome}/bin/larmord_extractor -rna -ring -cutoffRing 10.0 -identification ${id} -csfile ${data} -beta -3 -cutoff $2 ${PDB}/model_${m}_minimize.pdb	
        [[ $m != "1" ]] && ${larmorDHome}/bin/larmord_extractor -rna -ring -cutoffRing 10.0 -identification ${id} -csfile ${data} -beta -3 -cutoff $2 ${PDB}/model_${m}_minimize.pdb  | grep -v frame

	done
fi
