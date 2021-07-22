#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex9'
testname='runex9_1'
label='mat_tutorials-ex9_1'
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
args=''
diff_args=''
timeoutfactor=1

mpiexec=${PETSCMPIEXEC:-"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/bin/mpiexec"}
diffexec=${PETSCDIFF:-"${petsc_bindir}/petscdiff"}

. "${config_dir}/petsc_harness.sh"

# The diff flags come from script arguments
diff_exe="${diffexec} ${diff_flags} ${diff_args}"
mpiexec="${mpiexec} ${mpiexec_flags}"





nsize_in=${nsize:-2}
mat_composite_merge_in=${mat_composite_merge:-0 1}
mat_composite_merge_mvctx_in=${mat_composite_merge_mvctx:-0 1}


for insize in ${nsize_in}; do
   for imat_composite_merge in ${mat_composite_merge_in}; do
      for imat_composite_merge_mvctx in ${mat_composite_merge_mvctx_in}; do

         petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args}  -mat_composite_merge ${imat_composite_merge} -mat_composite_merge_mvctx ${imat_composite_merge_mvctx}" ex9_1.tmp ${testname}.err "${label}+mat_composite_merge-${imat_composite_merge}_mat_composite_merge_mvctx-${imat_composite_merge_mvctx}" 
         res=$?

         if test $res = 0; then
            petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/mat/tutorials/output/ex9_1.out ex9_1.tmp" diff-${testname}.out diff-${testname}.out diff-${label}+mat_composite_merge-${imat_composite_merge}_mat_composite_merge_mvctx-${imat_composite_merge_mvctx} ""
         else
            petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
         fi

      done
   done
done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
