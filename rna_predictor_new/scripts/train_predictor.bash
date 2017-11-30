#!/bin/bash
if [[ $# -ne 3 ]]
then
	echo "Script to calculate chemical shifts"
	echo "usage $0: <nucleus> <struct> <r_cut>"
else
	source ~/.bashrc	
	
  nucleus=$1
  struct=$2
  r_cut=$3
  
  # script to train chemical shift predictor for a given nucleus
  awk -v nucleus=${nucleus} '{if($1==nucleus) print}' data/initial_parameters.txt > data/initial_parameters_tmp.txt
  cat ????/train_features_${struct}_${r_cut}.txt | awk -v nucleus=${nucleus} '{ if($5==nucleus) print}' > data/train_tmp.txt   
  ./scripts/train.R --verbose --ga_populations=68 --ga_cycles=10000 --output_parameters=data/ga_model_${nucleus}.txt data/train_tmp.txt
fi
