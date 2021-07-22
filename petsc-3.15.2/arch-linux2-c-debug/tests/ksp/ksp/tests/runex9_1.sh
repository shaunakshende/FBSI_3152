#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex9'
testname='runex9_1'
label='ksp_ksp_tests-ex9_1'
runfiles=''
wPETSC_DIR='/home/shaunak/Desktop/petsc-3.15.2'
petsc_dir='/home/shaunak/Desktop/petsc-3.15.2'
petsc_arch='arch-linux2-c-debug'
# Must be consistent with gmakefile.test
testlogtapfile=/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests/test_${petsc_arch}_tap.log
testlogerrfile=/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests/test_${petsc_arch}_err.log
config_dir='/home/shaunak/Desktop/petsc-3.15.2/config'
filter='sed -e "s/CONVERGED_RTOL/CONVERGED_ATOL/g"'
filter_output=''
petsc_bindir='/home/shaunak/Desktop/petsc-3.15.2/lib/petsc/bin'
DATAFILESPATH=${DATAFILESPATH:-""}
args='-ksp_converged_reason -ksp_error_if_not_converged'
diff_args=''
timeoutfactor=1

mpiexec=${PETSCMPIEXEC:-"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/bin/mpiexec"}
diffexec=${PETSCDIFF:-"${petsc_bindir}/petscdiff"}

. "${config_dir}/petsc_harness.sh"

# The diff flags come from script arguments
diff_exe="${diffexec} ${diff_flags} ${diff_args}"
mpiexec="${mpiexec} ${mpiexec_flags}"





nsize_in=${nsize:-1}
pc_fieldsplit_diag_use_amat_in=${pc_fieldsplit_diag_use_amat:-0 1}
pc_fieldsplit_diag_use_amat_in=${pc_fieldsplit_diag_use_amat:-0 1}
pc_fieldsplit_type_in=${pc_fieldsplit_type:-additive multiplicative}


for insize in ${nsize_in}; do
   for ipc_fieldsplit_diag_use_amat in ${pc_fieldsplit_diag_use_amat_in}; do
      for ipc_fieldsplit_diag_use_amat in ${pc_fieldsplit_diag_use_amat_in}; do
         for ipc_fieldsplit_type in ${pc_fieldsplit_type_in}; do

            petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args}  -pc_fieldsplit_diag_use_amat ${ipc_fieldsplit_diag_use_amat} -pc_fieldsplit_diag_use_amat ${ipc_fieldsplit_diag_use_amat} -pc_fieldsplit_type ${ipc_fieldsplit_type}" ex9_1.tmp ${testname}.err "${label}+pc_fieldsplit_diag_use_amat-${ipc_fieldsplit_diag_use_amat}_pc_fieldsplit_diag_use_amat-${ipc_fieldsplit_diag_use_amat}_pc_fieldsplit_type-${ipc_fieldsplit_type}" 
            res=$?

            if test $res = 0; then
               petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/ksp/ksp/tests/output/ex9_1.out ex9_1.tmp" diff-${testname}.out diff-${testname}.out diff-${label}+pc_fieldsplit_diag_use_amat-${ipc_fieldsplit_diag_use_amat}_pc_fieldsplit_diag_use_amat-${ipc_fieldsplit_diag_use_amat}_pc_fieldsplit_type-${ipc_fieldsplit_type} ""
            else
               petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
            fi

         done
      done
   done
done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
