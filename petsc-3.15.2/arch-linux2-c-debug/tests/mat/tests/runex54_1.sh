#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex54'
testname='runex54_1'
label='mat_tests-ex54_1'
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





nsize_in=${nsize:-1 3}
mat_block_size_in=${mat_block_size:-1 3 4 6 8}
ov_in=${ov:-1 3}
mat_size_in=${mat_size:-11 13}
nd_in=${nd:-7}


for insize in ${nsize_in}; do
   for imat_block_size in ${mat_block_size_in}; do
      for iov in ${ov_in}; do
         for imat_size in ${mat_size_in}; do
            for ind in ${nd_in}; do

               petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args}  -mat_block_size ${imat_block_size} -ov ${iov} -mat_size ${imat_size} -nd ${ind}" ex54_1.tmp ${testname}.err "${label}+nsize-${insize}mat_block_size-${imat_block_size}_ov-${iov}_mat_size-${imat_size}_nd-${ind}" 
               res=$?

               if test $res = 0; then
                  petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/mat/tests/output/ex54.out ex54_1.tmp" diff-${testname}.out diff-${testname}.out diff-${label}+nsize-${insize}mat_block_size-${imat_block_size}_ov-${iov}_mat_size-${imat_size}_nd-${ind} ""
               else
                  petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
               fi

            done
         done
      done
   done
done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
