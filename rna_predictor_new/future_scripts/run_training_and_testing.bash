#!/bin/bash
if [[ $# -ne 1 ]]
then
	echo "Execute the training and testing of the new larmord parameters"
	echo "usage $0: <nucleus>" 
else
	source ~/.bashrc	
  nucleus=$1
  bash train_predictor.bash "${nucleus}" "data/train_features.txt" "data/larmord_features_train.txt"
  bash test_predictor.bash "${nucleus}" "data/larmord_features_test.txt" "data/bayesian_ga_model_${nucleus}.txt" "data/initial_parameters.txt"
fi