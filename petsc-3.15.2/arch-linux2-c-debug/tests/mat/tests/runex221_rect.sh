#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex221'
testname='runex221_rect'
label='mat_tests-ex221_rect'
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
args='-loop 3'
diff_args=''
timeoutfactor=1

mpiexec=${PETSCMPIEXEC:-"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/bin/mpiexec"}
diffexec=${PETSCDIFF:-"${petsc_bindir}/petscdiff"}

. "${config_dir}/petsc_harness.sh"

# The diff flags come from script arguments
diff_exe="${diffexec} ${diff_flags} ${diff_args}"
mpiexec="${mpiexec} ${mpiexec_flags}"





nsize_in=${nsize:-1 3}
keep_in=${keep:-0 1}
M_in=${M:-12 19}
N_in=${N:-19 12}
submat_in=${submat:-0 1}
test_axpy_different_in=${test_axpy_different:-0 1}


for insize in ${nsize_in}; do
   for ikeep in ${keep_in}; do
      for iM in ${M_in}; do
         for iN in ${N_in}; do
            for isubmat in ${submat_in}; do
               for itest_axpy_different in ${test_axpy_different_in}; do

                  petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args}  -keep ${ikeep} -M ${iM} -N ${iN} -submat ${isubmat} -test_axpy_different ${itest_axpy_different}" ex221_rect.tmp ${testname}.err "${label}+nsize-${insize}keep-${ikeep}_M-${iM}_N-${iN}_submat-${isubmat}_test_axpy_different-${itest_axpy_different}" 
                  res=$?

                  if test $res = 0; then
                     petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/mat/tests/output/ex221_1.out ex221_rect.tmp" diff-${testname}.out diff-${testname}.out diff-${label}+nsize-${insize}keep-${ikeep}_M-${iM}_N-${iN}_submat-${isubmat}_test_axpy_different-${itest_axpy_different} ""
                  else
                     petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
                  fi

               done
            done
         done
      done
   done
done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
