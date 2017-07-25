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
	data=${larmorDHome}/rna_predictor/data/measured_shifts_corrected_${PDB}.dat
	
	# goto working directory
	cd 	${larmorDHome}/rna_predictor/
	
	# Get number of models
	models=`ls ${PDB}/*pdb | wc -l | awk '{print $1}'`
	models=`seq 1 ${models}`	
  
  # Loop over NMR bundle and extract structure features
  for m in $models
  do
  	id=`echo "${PDB}_${m}"`
  	[[ $m == "1" ]] && ${larmorDHome}/bin/larmord_extractor -rna -ring -cutoffRing 10.0 -identification ${id} -csfile ${data} -beta -3 -cutoff 9999 ${PDB}/model_${m}.pdb
  	[[ $m != "1" ]] && ${larmorDHome}/bin/larmord_extractor -rna -ring -cutoffRing 10.0 -identification ${id} -csfile ${data} -beta -3 -cutoff 9999 ${PDB}/model_${m}.pdb  | grep -v frame	
  done
fi
