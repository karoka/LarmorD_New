#!/bin/bash
if [[ $# -ne 1 ]]
then
	echo "Script to calculate chemical shifts"
	echo "usage $0: <nucleus>"
else
	source ~/.bashrc	
	
  nucleus=$1
  
  # script to train chemical shift predictor for a given nucleus
  cat ????/train_features_min_15.0A.txt | awk -v nucleus=${nucleus} '{ if($5==nucleus) print}' > data/train_tmp.txt   
  ./scripts/train.R --verbose --ga_populations=68 --ga_cycles=1000 --output_parameters=data/ga_model_${nucleus}.txt data/train_tmp.txt
fi