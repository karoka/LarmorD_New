#!/bin/bash

LarmorD_home=/Users/blairwhittington/Documents/LarmorD_New/rna_predictor

# Extract features for training set RNAs
i=1
cut_offs="5.0 10.0 15.0 20.0 25.0"
rnas="1L1W 1LC6 1LDZ 1NC0 1OW9 1PJY 1R7W 1R7Z 1UUU 1YSV 2FDT 2GM0 2JXQ 2JXS 2K3Z 2K41 2KOC 2KYD 2LBJ 2LBL 2LDL 2LDT 2LHP 2LI4 2LK3 2LP9 2LPA 2LU0 2LUB 2LV0 2RN1 2Y95 4A4S 4A4T 4A4U"
#rnas="1KKA 1L1W 1LC6 1LDZ 1NC0 1OW9 1PJY 1R7W 1R7Z 1UUU 1YSV 2FDT 2GM0 2JXQ 2JXS 2K3Z 2K41 2KOC 2KYD 2LBJ 2LBL 2LDL 2LDT 2LHP 2LI4 2LK3 2LP9 2LPA 2LU0 2LUB 2LV0 2RN1 2Y95 4A4S 4A4T 4A4U"

for cut_off in $cut_offs
do
#rm -f  train_features_unmin_${cut_off}A.txt 
rm -f  train_features_min_${cut_off}A.txt
	for rna in $rnas
	do
		#[[ $i == 1 ]] && bash ${LarmorD_home}/extract_features.bash $rna $cut_off | tee ${LarmorD_home}/${rna}/train_features_unmin_${cut_off}A.txt 
		#[[ $i != 1 ]] && bash ${LarmorD_home}/extract_features.bash $rna $cut_off | grep -v frame | tee -a ${LarmorD_home}/${rna}/train_features_unmin_${cut_off}A.txt	
		#i=999	
	
		[[ $i == 1 ]] && bash ${LarmorD_home}/extract_features.bash $rna $cut_off | tee ${LarmorD_home}/${rna}/train_features_min_${cut_off}A.txt
		[[ $i != 1 ]] && bash ${LarmorD_home}/extract_features.bash $rna $cut_off | grep -v frame | tee -a ${LarmorD_home}/${rna}/train_features_min_${cut_off}A.txt
		i=999
	done
done
 
 
# Extract features for testing set RNAs
i=1
cut_offs="5.0 10.0 15.0 20.0 25.0"
rnas="2JWV 2K63 2K64 2K65 2K66 2LPS 2LQZ 2LUN 2LX1 2M12 2M21 2M22 2M23 2M24 2M57 2M8K 2MEQ 2MHI 2MI0 2MIS 1SCL 1XHP 1ZC5 2QH2 2QH3 2QH4"
#rnas="1Z2J 2JWV 2K63 2K64 2K65 2K66 2LPS 2LQZ 2LUN 2LX1 2M12 2M21 2M22 2M23 2M24 2M57 2M8K 2MEQ 2MHI 2MI0 2MIS 1SCL 1XHP 1ZC5 2QH2 2QH3 2QH4"

for cut_off in $cut_offs
do
#rm -f  test_features_unmin_${cut_off}A.txt
rm -f test_features_min_${cut_off}A.txt
	for rna in $rnas
	do
		#[[ $i == 1 ]] && bash ${LarmorD_home}/extract_features.bash $rna $cut_off | tee ${LarmorD_home}/${rna}/test_features_unmin_${cut_off}A.txt 
		#[[ $i != 1 ]] && bash ${LarmorD_home}/extract_features.bash $rna $cut_off | grep -v frame | tee -a ${LarmorD_home}/${rna}/test_features_unmin_${cut_off}A.txt
		#i=999

		[[ $i == 1 ]] && bash ${LarmorD_home}/extract_features.bash $rna $cut_off | tee ${LarmorD_home}/${rna}/test_features_min_${cut_off}A.txt
		[[ $i != 1 ]] && bash ${LarmorD_home}/extract_features.bash $rna $cut_off | grep -v frame | tee -a ${LarmorD_home}/${rna}/test_features_min_${cut_off}A.txt
		i=999
	done
done

#awk '{ if(NF==80) print}' test_features.txt > tmp.txt; mv tmp.txt test_features.txt
#awk '{ if(NF==80) print}' train_features.txt > tmp.txt; mv tmp.txt train_features.txt

