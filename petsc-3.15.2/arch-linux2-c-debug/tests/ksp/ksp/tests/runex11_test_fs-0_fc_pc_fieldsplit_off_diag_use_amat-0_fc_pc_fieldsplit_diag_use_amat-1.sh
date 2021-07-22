#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex11'
testname='runex11_test_fs-0_fc_pc_fieldsplit_off_diag_use_amat-0_fc_pc_fieldsplit_diag_use_amat-1'
label='ksp_ksp_tests-ex11_test_fs-0_fc_pc_fieldsplit_off_diag_use_amat-0_fc_pc_fieldsplit_diag_use_amat-1'
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
args='-f ${DATAFILESPATH}/matrices/underworld32.gz -fc_ksp_view -fc_ksp_monitor_short -fc_ksp_type fgmres -fc_ksp_max_it 4000 -fc_ksp_diagonal_scale -fc_pc_type fieldsplit -fc_pc_fieldsplit_type SCHUR -fc_pc_fieldsplit_schur_fact_type UPPER -fc_fieldsplit_velocity_ksp_type cg -fc_fieldsplit_velocity_pc_type cholesky -fc_fieldsplit_velocity_pc_factor_mat_ordering_type nd -fc_fieldsplit_pressure_ksp_max_it 100 -fc_fieldsplit_pressure_ksp_constant_null_space -fc_fieldsplit_pressure_ksp_monitor_short -fc_fieldsplit_pressure_pc_type lsc -fc_fieldsplit_pressure_lsc_ksp_type cg -fc_fieldsplit_pressure_lsc_ksp_max_it 100 -fc_fieldsplit_pressure_lsc_ksp_constant_null_space -fc_fieldsplit_pressure_lsc_ksp_converged_reason -fc_fieldsplit_pressure_lsc_pc_type icc -test_fs 0 -fc_pc_fieldsplit_off_diag_use_amat 0 -fc_pc_fieldsplit_diag_use_amat 1'
diff_args=''
timeoutfactor=1

mpiexec=${PETSCMPIEXEC:-"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/bin/mpiexec"}
diffexec=${PETSCDIFF:-"${petsc_bindir}/petscdiff"}

. "${config_dir}/petsc_harness.sh"

# The diff flags come from script arguments
diff_exe="${diffexec} ${diff_flags} ${diff_args}"
mpiexec="${mpiexec} ${mpiexec_flags}"



if ! $force; then
    if test -z "${DATAFILESPATH}"; then
        petsc_report_tapoutput "" "${label}" "SKIP Requires DATAFILESPATH"
        total=1; skip=1
        petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
        exit
    fi
fi


nsize_in=${nsize:-1}


for insize in ${nsize_in}; do

   petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args} " ex11_test_fs-0_fc_pc_fieldsplit_off_diag_use_amat-0_fc_pc_fieldsplit_diag_use_amat-1.tmp ${testname}.err "${label}" 
   res=$?

   if test $res = 0; then
      petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/ksp/ksp/tests/output/ex11_test_fs-0_fc_pc_fieldsplit_off_diag_use_amat-0_fc_pc_fieldsplit_diag_use_amat-1.out ex11_test_fs-0_fc_pc_fieldsplit_off_diag_use_amat-0_fc_pc_fieldsplit_diag_use_amat-1.tmp" diff-${testname}.out diff-${testname}.out diff-${label} ""
   else
      petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
   fi

done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
