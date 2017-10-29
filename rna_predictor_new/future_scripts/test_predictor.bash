#!/bin/bash
if [[ $# -ne 4 ]]
then
	echo "Script to evaluate the performance of the refined larmord parameters"
	echo "usage $0: <nucleus> <path/to/decoy/data/for/genetic/algorithm/optimization> <path/to/bayesian/ga/parameters> <path/to/orginal/larmord/parameters>"
else
	source ~/.bashrc	
  nucleus=$1
  decoy_data=$2
  new_parameter_data=$3
  larmord_parameter_data=$4
  
  # filter data
  awk -v nucleus=${nucleus} '{ if($7==nucleus) print}' ${decoy_data} > data/decoy_tmp.txt 
  awk -v nucleus=${nucleus} '{ if($1==nucleus) print}' ${larmord_parameter_data} > data/larmord_tmp.txt 
  
  # script to train chemical shift predictor for a given nucleus
  ./test_bayes_ga.R --comparison_fparameters=data/larmord_tmp.txt ${new_parameter_data} data/decoy_tmp.txt 
fi