cd tests
../bin/larmord  -csfile measured_shifts_A003.dat \
-parmfile ../data/larmorD_alphas_betas.dat \
-reffile ../data/larmorD_reference_shifts.dat \
-accfile ../data/larmorD_accuracy_weight.dat \
-cutoff 15.0 A003.pdb
