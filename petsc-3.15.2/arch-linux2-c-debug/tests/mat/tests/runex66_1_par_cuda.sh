#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex66'
testname='runex66_1_par_cuda'
label='mat_tests-ex66_1_par_cuda'
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
args='-n 32 -kernel 1 -dim 1 -ldc 12'
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


nsize_in=${nsize:-2}
testlayout_in=${testlayout:-0 1}
bgpu_in=${bgpu:-0 1}
cgpu_in=${cgpu:-0 1}


for insize in ${nsize_in}; do
   for itestlayout in ${testlayout_in}; do
      for ibgpu in ${bgpu_in}; do
         for icgpu in ${cgpu_in}; do

            petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args}  -testlayout ${itestlayout} -bgpu ${ibgpu} -cgpu ${icgpu}" ex66_1_par_cuda.tmp ${testname}.err "${label}+testlayout-${itestlayout}_bgpu-${ibgpu}_cgpu-${icgpu}" 
            res=$?

            if test $res = 0; then
               petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/mat/tests/output/ex66_1_par.out ex66_1_par_cuda.tmp" diff-${testname}.out diff-${testname}.out diff-${label}+testlayout-${itestlayout}_bgpu-${ibgpu}_cgpu-${icgpu} ""
            else
               petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
            fi

         done
      done
   done
done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
