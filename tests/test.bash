cd tests
../bin/larmord  -csfile measured_shifts_A003.dat \
-parmfile ../data/larmorD_alphas_betas.dat \
-reffile ../data/larmorD_reference_shifts.dat \
-accfile ../data/larmorD_accuracy_weight.dat \
-cutoff 15.0 A003.pdb
echo " "
../bin/larmord  -csfile measured_shifts_2ADT.dat \
-parmfile ../data/larmorD_alphas_betas_rna.dat \
-reffile ../data/larmorD_reference_shifts_rna.dat \
-accfile ../data/larmorD_accuracy_weight_rna.dat \
-cutoff 9999.0 2ADT.pdb
echo " "
../bin/larmord  -csfile measured_shifts_2ADT.dat \
-parmfile ../data/larmorD_alphas_betas_rna.dat \
-reffile ../data/larmorD_reference_shifts_rna.dat \
-accfile ../data/larmorD_accuracy_weight_rna.dat \
-printError \
-cutoff 9999.0 2ADT.pdb

