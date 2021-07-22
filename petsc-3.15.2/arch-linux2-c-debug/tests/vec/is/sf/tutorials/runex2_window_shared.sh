#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex2'
testname='runex2_window_shared'
label='vec_is_sf_tutorials-ex2_window_shared'
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
args='-sf_type window -sf_window_flavor shared'
diff_args=''
timeoutfactor=1

mpiexec=${PETSCMPIEXEC:-"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/bin/mpiexec"}
diffexec=${PETSCDIFF:-"${petsc_bindir}/petscdiff"}

. "${config_dir}/petsc_harness.sh"

# The diff flags come from script arguments
diff_exe="${diffexec} ${diff_flags} ${diff_args}"
mpiexec="${mpiexec} ${mpiexec_flags}"



if ! $force; then
    petsc_report_tapoutput "" "${label}" "SKIP Null requirement not met: define(PETSC_HAVE_MPICH_NUMVERSION)"
    total=1; skip=1
    petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
    exit
fi


nsize_in=${nsize:-2}
sf_window_sync_in=${sf_window_sync:-fence active lock}


for insize in ${nsize_in}; do
   for isf_window_sync in ${sf_window_sync_in}; do

      petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args}  -sf_window_sync ${isf_window_sync}" ex2_window_shared.tmp ${testname}.err "${label}+sf_window_sync-${isf_window_sync}" 
      res=$?

      if test $res = 0; then
         petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/vec/is/sf/tutorials/output/ex2_basic.out ex2_window_shared.tmp" diff-${testname}.out diff-${testname}.out diff-${label}+sf_window_sync-${isf_window_sync} ""
      else
         petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
      fi

   done
done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
