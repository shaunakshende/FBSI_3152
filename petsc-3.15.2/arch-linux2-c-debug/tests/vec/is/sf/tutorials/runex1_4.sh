#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex1'
testname='runex1_4'
label='vec_is_sf_tutorials-ex1_4'
runfiles=''
wPETSC_DIR='/home/shaunak/Desktop/petsc-3.15.2'
petsc_dir='/home/shaunak/Desktop/petsc-3.15.2'
petsc_arch='arch-linux2-c-debug'
# Must be consistent with gmakefile.test
testlogtapfile=/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests/test_${petsc_arch}_tap.log
testlogerrfile=/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests/test_${petsc_arch}_err.log
config_dir='/home/shaunak/Desktop/petsc-3.15.2/config'
filter='grep -v "type" | grep -v "sort"'
filter_output=''
petsc_bindir='/home/shaunak/Desktop/petsc-3.15.2/lib/petsc/bin'
DATAFILESPATH=${DATAFILESPATH:-""}
args='-test_gather -sf_type window'
diff_args=''
timeoutfactor=1

mpiexec=${PETSCMPIEXEC:-"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/bin/mpiexec"}
diffexec=${PETSCDIFF:-"${petsc_bindir}/petscdiff"}

. "${config_dir}/petsc_harness.sh"

# The diff flags come from script arguments
diff_exe="${diffexec} ${diff_flags} ${diff_args}"
mpiexec="${mpiexec} ${mpiexec_flags}"





nsize_in=${nsize:-4}
sf_window_sync_in=${sf_window_sync:-fence active lock}
sf_window_flavor_in=${sf_window_flavor:-create dynamic allocate}


for insize in ${nsize_in}; do
   for isf_window_sync in ${sf_window_sync_in}; do
      for isf_window_flavor in ${sf_window_flavor_in}; do

         petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args}  -sf_window_sync ${isf_window_sync} -sf_window_flavor ${isf_window_flavor}" ex1_4.tmp ${testname}.err "${label}+sf_window_sync-${isf_window_sync}_sf_window_flavor-${isf_window_flavor}" 
         res=$?

         if test $res = 0; then
            petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/vec/is/sf/tutorials/output/ex1_4.out ex1_4.tmp" diff-${testname}.out diff-${testname}.out diff-${label}+sf_window_sync-${isf_window_sync}_sf_window_flavor-${isf_window_flavor} ""
         else
            petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
         fi

      done
   done
done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
