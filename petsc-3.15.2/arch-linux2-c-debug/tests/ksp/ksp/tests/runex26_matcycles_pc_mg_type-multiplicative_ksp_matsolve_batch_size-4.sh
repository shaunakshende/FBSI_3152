#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex26'
testname='runex26_matcycles_pc_mg_type-multiplicative_ksp_matsolve_batch_size-4'
label='ksp_ksp_tests-ex26_matcycles_pc_mg_type-multiplicative_ksp_matsolve_batch_size-4'
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
args='-ksp_view_final_residual -ksp_type preonly -pc_type mg -mx 5 -my 5 -pc_mg_levels 3 -pc_mg_galerkin -ksp_monitor -mg_levels_ksp_type richardson -mg_levels_pc_type jacobi -nrhs 7 -pc_mg_type multiplicative -ksp_matsolve_batch_size 4'
diff_args=''
timeoutfactor=1

mpiexec=${PETSCMPIEXEC:-"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/bin/mpiexec"}
diffexec=${PETSCDIFF:-"${petsc_bindir}/petscdiff"}

. "${config_dir}/petsc_harness.sh"

# The diff flags come from script arguments
diff_exe="${diffexec} ${diff_flags} ${diff_args}"
mpiexec="${mpiexec} ${mpiexec_flags}"





nsize_in=${nsize:-1 2}


for insize in ${nsize_in}; do

   petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args} " ex26_matcycles_pc_mg_type-multiplicative_ksp_matsolve_batch_size-4.tmp ${testname}.err "${label}+nsize-${insize}" 
   res=$?

   if test $res = 0; then
      petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/ksp/ksp/tests/output/ex26_matcycles_pc_mg_type-multiplicative_ksp_matsolve_batch_size-4.out ex26_matcycles_pc_mg_type-multiplicative_ksp_matsolve_batch_size-4.tmp" diff-${testname}.out diff-${testname}.out diff-${label}+nsize-${insize} ""
   else
      petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
   fi

done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
