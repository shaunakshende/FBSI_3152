#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex27'
testname='runex27_6_hdf5_dense'
label='ksp_ksp_tutorials-ex27_6_hdf5_dense'
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
args='-ksp_converged_reason -ksp_monitor_short -ksp_rtol 1e-5 -ksp_max_it 10 -solve_normal 0 -ksp_type lsqr -hdf5 -x0_name x -f ${DATAFILESPATH}/matrices/matlab/small_dense.mat -mat_type dense'
diff_args=''
timeoutfactor=1

mpiexec=${PETSCMPIEXEC:-"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/bin/mpiexec"}
diffexec=${PETSCDIFF:-"${petsc_bindir}/petscdiff"}

. "${config_dir}/petsc_harness.sh"

# The diff flags come from script arguments
diff_exe="${diffexec} ${diff_flags} ${diff_args}"
mpiexec="${mpiexec} ${mpiexec_flags}"



if ! $force; then
    petsc_report_tapoutput "" "${label}" "SKIP Requires DATAFILESPATH, PETSC_HAVE_HDF5 requirement not met, Required: define(PETSC_HDF5_HAVE_ZLIB)"
    total=1; skip=1
    petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
    exit
fi


nsize_in=${nsize:-1 2 4 8}
test_custom_layout_in=${test_custom_layout:-0 1}


for insize in ${nsize_in}; do
   for itest_custom_layout in ${test_custom_layout_in}; do

      petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args}  -test_custom_layout ${itest_custom_layout}" ex27_6_hdf5_dense.tmp ${testname}.err "${label}+nsize-${insize}test_custom_layout-${itest_custom_layout}" 
      res=$?

      if test $res = 0; then
         petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/ksp/ksp/tutorials/output/ex27_6_hdf5_dense.out ex27_6_hdf5_dense.tmp" diff-${testname}.out diff-${testname}.out diff-${label}+nsize-${insize}test_custom_layout-${itest_custom_layout} ""
      else
         petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
      fi

   done
done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
