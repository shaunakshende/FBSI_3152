#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex66'
testname='runex66_2_cuda'
label='mat_tests-ex66_2_cuda'
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
args='-kernel 0 -dim 2 -symm 1 -checkexpl'
diff_args=''
timeoutfactor=1

mpiexec=${PETSCMPIEXEC:-"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/bin/mpiexec"}
diffexec=${PETSCDIFF:-"${petsc_bindir}/petscdiff"}

. "${config_dir}/petsc_harness.sh"

# The diff flags come from script arguments
diff_exe="${diffexec} ${diff_flags} ${diff_args}"
mpiexec="${mpiexec} ${mpiexec_flags}"



if ! $force; then
    petsc_report_tapoutput "" "${label}" "SKIP PETSC_HAVE_HARA requirement not met, PETSC_HAVE_HARA requirement not met, PETSC_HAVE_CUDA requirement not met, PETSC_HAVE_HARA requirement not met"
    total=1; skip=1
    petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
    exit
fi


nsize_in=${nsize:-1}
n_in=${n:-17 33}
bgpu_in=${bgpu:-0 1}
agpu_in=${agpu:-0 1}


for insize in ${nsize_in}; do
   for in in ${n_in}; do
      for ibgpu in ${bgpu_in}; do
         for iagpu in ${agpu_in}; do

            petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args}  -n ${in} -bgpu ${ibgpu} -agpu ${iagpu}" ex66_2_cuda.tmp ${testname}.err "${label}+n-${in}_bgpu-${ibgpu}_agpu-${iagpu}" 
            res=$?

            if test $res = 0; then
               petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/mat/tests/output/ex66_2.out ex66_2_cuda.tmp" diff-${testname}.out diff-${testname}.out diff-${label}+n-${in}_bgpu-${ibgpu}_agpu-${iagpu} ""
            else
               petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
            fi

         done
      done
   done
done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
