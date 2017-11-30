#/bin/bash

home=/home/blairjw/franklab/new_larmord_2.0/LarmorD_New/rna_predictor_new

rnas="1KKA 1L1W 1LC6 1LDZ 1NC0 1OW9 1PJY 1R7W 1R7Z 1UUU 1YSV 2FDT 2GM0 2JXS 2K3Z 2K41 2KOC 2KYD 2LBJ 2LBL 2LDL 2LDT 2LHP 2LI4 2LK3 2LP9 2LPA 2LU0 2LUB 2LV0 2RN1 2Y95 4A4S 4A4T 4A4U" 

i=0
for rna in $rnas
do
	i=$(($i+1))
	echo "processing $i: $rna"
	cd ${home}/${rna}

	echo "#!/bin/bash">job_temp.sh
        echo "#SBATCH -A frank" >>job_temp.sh
        echo "#SBATCH --job-name=${rna}" >>job_temp.sh
        echo "#SBATCH --time=3:00:00" >>job_temp.sh
        echo "#SBATCH --ntasks=1" >>job_temp.sh
        echo "#SBATCH --nodes=1" >>job_temp.sh
        echo "" >>job_temp.sh
        echo "" >>job_temp.sh
        #echo "module load rosetta" >>job_temp.sh
        #echo "export ROSETTA=/export/apps/rosetta/src/rosetta_src_2016.32.58837_bundle" >>job_temp.sh
        echo "source ~/.bashrc" >>job_temp.sh
	
	sbatch -A frank -p frank job_temp.sh

	cd ${home}
done 
