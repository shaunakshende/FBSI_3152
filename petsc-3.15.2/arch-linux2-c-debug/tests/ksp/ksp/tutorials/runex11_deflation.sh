#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex11'
testname='runex11_deflation'
label='ksp_ksp_tutorials-ex11_deflation'
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
args='-norandom -pc_type deflation -ksp_monitor_short'
diff_args=''
timeoutfactor=1

mpiexec=${PETSCMPIEXEC:-"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/bin/mpiexec"}
diffexec=${PETSCDIFF:-"${petsc_bindir}/petscdiff"}

. "${config_dir}/petsc_harness.sh"

# The diff flags come from script arguments
diff_exe="${diffexec} ${diff_flags} ${diff_args}"
mpiexec="${mpiexec} ${mpiexec_flags}"



if ! $force; then
    petsc_report_tapoutput "" "${label}" "SKIP PETSC_USE_COMPLEX requirement not met, PETSC_HAVE_SUPERLU_DIST requirement not met, PETSC_USE_COMPLEX requirement not met"
    total=1; skip=1
    petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
    exit
fi


nsize_in=${nsize:-6}


for insize in ${nsize_in}; do

   petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args}  -norandom -pc_type deflation -ksp_monitor_short" ex11_deflation.tmp ${testname}.err "${label}+a" 
   res=$?

   if test $res = 0; then
      petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/ksp/ksp/tutorials/output/ex11_deflation.out ex11_deflation.tmp" diff-${testname}.out diff-${testname}.out diff-${label}+a ""
   else
      petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
   fi

done

nsize_in=${nsize:-3}
pc_deflation_compute_space_in=${pc_deflation_compute_space:-db2 aggregation}


for insize in ${nsize_in}; do
   for ipc_deflation_compute_space in ${pc_deflation_compute_space_in}; do

      petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args}  -pc_deflation_compute_space ${ipc_deflation_compute_space} " ex11_deflation.tmp ${testname}.err "${label}+b+pc_deflation_compute_space-${ipc_deflation_compute_space}" 
      res=$?

      if test $res = 0; then
         petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/ksp/ksp/tutorials/output/ex11_deflation.out ex11_deflation.tmp" diff-${testname}.out diff-${testname}.out diff-${label}+b+pc_deflation_compute_space-${ipc_deflation_compute_space} ""
      else
         petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
      fi

   done
done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
