#!/bin/bash
#SBATCH -A frank
#SBATCH --job-name=1PJY
#SBATCH --time=3:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1


module load rosetta
export ROSETTA=/export/apps/rosetta/src/rosetta_src_2016.32.58837_bundle
source ~/.bashrc
bash /export/nVerde/users/blairjw/frank_lab/larmord_2.0/LarmordD_New/minimize_for_predictor.bash 1PJY
