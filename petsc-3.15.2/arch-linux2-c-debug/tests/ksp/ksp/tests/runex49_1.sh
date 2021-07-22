#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex49'
testname='runex49_1'
label='ksp_ksp_tests-ex49_1'
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
args='-pc_type cholesky -herm 0'
diff_args=''
timeoutfactor=1

mpiexec=${PETSCMPIEXEC:-"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/bin/mpiexec"}
diffexec=${PETSCDIFF:-"${petsc_bindir}/petscdiff"}

. "${config_dir}/petsc_harness.sh"

# The diff flags come from script arguments
diff_exe="${diffexec} ${diff_flags} ${diff_args}"
mpiexec="${mpiexec} ${mpiexec_flags}"





nsize_in=${nsize:-1}
mat_type_in=${mat_type:-aij baij sbaij}
bs_in=${bs:-1 2 3 4 5 6 7 8 9 10 11 12}
conv_in=${conv:-0 1}


for insize in ${nsize_in}; do
   for imat_type in ${mat_type_in}; do
      for ibs in ${bs_in}; do
         for iconv in ${conv_in}; do

            petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args}  -mat_type ${imat_type} -bs ${ibs} -conv ${iconv}" ex49_1.tmp ${testname}.err "${label}+mat_type-${imat_type}_bs-${ibs}_conv-${iconv}" 
            res=$?

            if test $res = 0; then
               petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/ksp/ksp/tests/output/ex49_1.out ex49_1.tmp > diff-${testname}-0.out 2> diff-${testname}-0.out || ${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/ksp/ksp/tests/output/ex49_1_alt.out ex49_1.tmp" diff-${testname}.out diff-${testname}.out diff-${label}+mat_type-${imat_type}_bs-${ibs}_conv-${iconv} ""
            else
               petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
            fi

         done
      done
   done
done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
