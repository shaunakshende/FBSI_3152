#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex81'
testname='runex81_w_explicit_schur_outer_fieldsplit_sizes-2__3'
label='ksp_ksp_tutorials-ex81_w_explicit_schur_outer_fieldsplit_sizes-2__3'
runfiles=''
wPETSC_DIR='/home/shaunak/Desktop/petsc-3.15.2'
petsc_dir='/home/shaunak/Desktop/petsc-3.15.2'
petsc_arch='arch-linux2-c-debug'
# Must be consistent with gmakefile.test
testlogtapfile=/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests/test_${petsc_arch}_tap.log
testlogerrfile=/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests/test_${petsc_arch}_err.log
config_dir='/home/shaunak/Desktop/petsc-3.15.2/config'
filter='sed -e "s/seq/mpi/g" -e "s/CONVERGED_ATOL/CONVERGED_RTOL/g" -e "s/ 1 MPI processes/ 3 MPI processes/g" -e "s/iterations [1-2]/iterations 1/g"'
filter_output=''
petsc_bindir='/home/shaunak/Desktop/petsc-3.15.2/lib/petsc/bin'
DATAFILESPATH=${DATAFILESPATH:-""}
args='-ksp_view_mat -pc_type fieldsplit -ksp_converged_reason -fieldsplit_pc_type jacobi -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type full -pc_fieldsplit_schur_precondition full -outer_fieldsplit_sizes 2,3'
diff_args=''
timeoutfactor=1

mpiexec=${PETSCMPIEXEC:-"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/bin/mpiexec"}
diffexec=${PETSCDIFF:-"${petsc_bindir}/petscdiff"}

. "${config_dir}/petsc_harness.sh"

# The diff flags come from script arguments
diff_exe="${diffexec} ${diff_flags} ${diff_args}"
mpiexec="${mpiexec} ${mpiexec_flags}"





nsize_in=${nsize:-1 3}


for insize in ${nsize_in}; do

   petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args} " ex81_w_explicit_schur_outer_fieldsplit_sizes-2__3.tmp ${testname}.err "${label}+nsize-${insize}" 
   res=$?

   if test $res = 0; then
      petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/ksp/ksp/tutorials/output/ex81_w_explicit_schur_outer_fieldsplit_sizes-2__3.out ex81_w_explicit_schur_outer_fieldsplit_sizes-2__3.tmp" diff-${testname}.out diff-${testname}.out diff-${label}+nsize-${insize} ""
   else
      petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
   fi

done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
