#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex3'
testname='runex3_2'
label='dm_impls_stag_tutorials-ex3_2'
runfiles=''
wPETSC_DIR='/home/shaunak/Desktop/petsc-3.15.2'
petsc_dir='/home/shaunak/Desktop/petsc-3.15.2'
petsc_arch='arch-linux2-c-debug'
# Must be consistent with gmakefile.test
testlogtapfile=/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests/test_${petsc_arch}_tap.log
testlogerrfile=/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests/test_${petsc_arch}_err.log
config_dir='/home/shaunak/Desktop/petsc-3.15.2/config'
filter=''
filter_output=''
petsc_bindir='/home/shaunak/Desktop/petsc-3.15.2/lib/petsc/bin'
DATAFILESPATH=${DATAFILESPATH:-""}
args='-ksp_monitor_short -ksp_converged_reason -pc_fieldsplit_schur_fact_type diag -fieldsplit_0_ksp_type preonly -fieldsplit_1_pc_type none -fieldsplit_0_pc_type gamg -fieldsplit_0_mg_levels_ksp_max_it 3 -fieldsplit_1_ksp_type gmres -fieldsplit_1_ksp_max_it 20 -fieldsplit_0_pc_gamg_esteig_ksp_type gmres'
diff_args=''
timeoutfactor=1

mpiexec=${PETSCMPIEXEC:-"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/bin/mpiexec"}
diffexec=${PETSCDIFF:-"${petsc_bindir}/petscdiff"}

. "${config_dir}/petsc_harness.sh"

# The diff flags come from script arguments
diff_exe="${diffexec} ${diff_flags} ${diff_args}"
mpiexec="${mpiexec} ${mpiexec_flags}"





nsize_in=${nsize:-4}


for insize in ${nsize_in}; do

   petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args} " ex3_2.tmp ${testname}.err "${label}" 
   res=$?

   if test $res = 0; then
      petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/dm/impls/stag/tutorials/output/ex3_2.out ex3_2.tmp > diff-${testname}-0.out 2> diff-${testname}-0.out || ${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/dm/impls/stag/tutorials/output/ex3_2_alt.out ex3_2.tmp > diff-${testname}-1.out 2> diff-${testname}-1.out || ${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/dm/impls/stag/tutorials/output/ex3_2_alt_2.out ex3_2.tmp > diff-${testname}-2.out 2> diff-${testname}-2.out || ${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/dm/impls/stag/tutorials/output/ex3_2_alt_3.out ex3_2.tmp > diff-${testname}-3.out 2> diff-${testname}-3.out || ${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/dm/impls/stag/tutorials/output/ex3_2_alt_4.out ex3_2.tmp" diff-${testname}.out diff-${testname}.out diff-${label} ""
   else
      petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
   fi

done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
