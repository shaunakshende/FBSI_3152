#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex31'
testname='runex31_0'
label='dm_impls_plex_tests-ex31_0'
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
args='-interpolate -use_initial_guess FALSE'
diff_args=''
timeoutfactor=1

mpiexec=${PETSCMPIEXEC:-"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/bin/mpiexec"}
diffexec=${PETSCDIFF:-"${petsc_bindir}/petscdiff"}

. "${config_dir}/petsc_harness.sh"

# The diff flags come from script arguments
diff_exe="${diffexec} ${diff_flags} ${diff_args}"
mpiexec="${mpiexec} ${mpiexec_flags}"



if ! $force; then
    petsc_report_tapoutput "" "${label}" "SKIP PETSC_HAVE_PARMETIS requirement not met"
    total=1; skip=1
    petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
    exit
fi


nsize_in=${nsize:-2 3 4}
faces_in=${faces:-2,3,4  5,4,3  7,11,5}
entity_depth_in=${entity_depth:-0 1}
parallel_in=${parallel:-FALSE TRUE}


for insize in ${nsize_in}; do
   for ifaces in ${faces_in}; do
      for ientity_depth in ${entity_depth_in}; do
         for iparallel in ${parallel_in}; do

            petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args}  -faces ${ifaces} -entity_depth ${ientity_depth} -parallel ${iparallel}" ex31_0.tmp ${testname}.err "${label}+nsize-${insize}faces-${ifaces}_entity_depth-${ientity_depth}_parallel-${iparallel}" 
            res=$?

            if test $res = 0; then
               petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/dm/impls/plex/tests/output/ex31_0.out ex31_0.tmp" diff-${testname}.out diff-${testname}.out diff-${label}+nsize-${insize}faces-${ifaces}_entity_depth-${ientity_depth}_parallel-${iparallel} ""
            else
               petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
            fi

         done
      done
   done
done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
