#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex70'
testname='runex70_2_cuda_real'
label='mat_tests-ex70_2_cuda_real'
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
args='-M 7 -N 9 -K 2 -testnest 0'
diff_args=''
timeoutfactor=1

mpiexec=${PETSCMPIEXEC:-"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/bin/mpiexec"}
diffexec=${PETSCDIFF:-"${petsc_bindir}/petscdiff"}

. "${config_dir}/petsc_harness.sh"

# The diff flags come from script arguments
diff_exe="${diffexec} ${diff_flags} ${diff_args}"
mpiexec="${mpiexec} ${mpiexec_flags}"



if ! $force; then
    petsc_report_tapoutput "" "${label}" "SKIP PETSC_HAVE_CUDA requirement not met"
    total=1; skip=1
    petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
    exit
fi


nsize_in=${nsize:-1}
local_in=${local:-0 1}
A_mat_type_in=${A_mat_type:-seqdensecuda seqdense}
xgpu_in=${xgpu:-0 1}
bgpu_in=${bgpu:-0 1}


for insize in ${nsize_in}; do
   for ilocal in ${local_in}; do
      for iA_mat_type in ${A_mat_type_in}; do
         for ixgpu in ${xgpu_in}; do
            for ibgpu in ${bgpu_in}; do

               petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args}  -local ${ilocal} -A_mat_type ${iA_mat_type} -xgpu ${ixgpu} -bgpu ${ibgpu}" ex70_2_cuda_real.tmp ${testname}.err "${label}+local-${ilocal}_A_mat_type-${iA_mat_type}_xgpu-${ixgpu}_bgpu-${ibgpu}" 
               res=$?

               if test $res = 0; then
                  petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/mat/tests/output/ex70_1.out ex70_2_cuda_real.tmp" diff-${testname}.out diff-${testname}.out diff-${label}+local-${ilocal}_A_mat_type-${iA_mat_type}_xgpu-${ixgpu}_bgpu-${ibgpu} ""
               else
                  petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
               fi

            done
         done
      done
   done
done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
