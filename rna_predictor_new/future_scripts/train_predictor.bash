#!/bin/bash
if [[ $# -ne 3 ]]
then
	echo "Script to calculate chemical shifts"
	echo "usage $0: <nucleus> <path/to/training/data/for/bayesian/regression> <path/to/decoy/data/for/genetic/algorithm/optimization>"
else
	source ~/.bashrc	
	
  nucleus=$1
  train_data=$2
  decoy_data=$3
  
  # script to train chemical shift predictor for a given nucleus
  awk -v nucleus=${nucleus} '{ if($5==nucleus) print}' ${train_data} > data/train_tmp.txt 
  awk -v nucleus=${nucleus} '{ if($7==nucleus) print}' ${decoy_data} > data/decoy_tmp.txt 
  
  ./train_bayes_ga.R --ga_cycles=250 --output_parameters=data/bayesian_ga_model_${nucleus}.txt data/train_tmp.txt data/decoy_tmp.txt 
fi