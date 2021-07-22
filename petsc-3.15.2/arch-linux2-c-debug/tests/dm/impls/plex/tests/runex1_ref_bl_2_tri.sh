#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex1'
testname='runex1_ref_bl_2_tri'
label='dm_impls_plex_tests-ex1_ref_bl_2_tri'
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
args='-dim 2 -domain_shape box -cell_simplex 1 -domain_box_sizes 5 -dm_view -interpolate -dm_plex_check_all 0 -dm_refine 1 -dm_plex_cell_refiner boundarylayer -ext_layers 3 -final_diagnostics -dm_plex_refine_boundarylayer_splits 4'
diff_args=''
timeoutfactor=1

mpiexec=${PETSCMPIEXEC:-"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/bin/mpiexec"}
diffexec=${PETSCDIFF:-"${petsc_bindir}/petscdiff"}

. "${config_dir}/petsc_harness.sh"

# The diff flags come from script arguments
diff_exe="${diffexec} ${diff_flags} ${diff_args}"
mpiexec="${mpiexec} ${mpiexec_flags}"



if ! $force; then
    petsc_report_tapoutput "" "${label}" "SKIP PETSC_HAVE_TRIANGLE requirement not met"
    total=1; skip=1
    petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
    exit
fi


nsize_in=${nsize:-1}
ext_hfirst_in=${ext_hfirst:-0 1}


for insize in ${nsize_in}; do
   for iext_hfirst in ${ext_hfirst_in}; do

      petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args}  -ext_hfirst ${iext_hfirst}" ex1_ref_bl_2_tri.tmp ${testname}.err "${label}+ext_hfirst-${iext_hfirst}" 
      res=$?

      if test $res = 0; then
         petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/dm/impls/plex/tests/output/ex1_ref_bl_2_tri.out ex1_ref_bl_2_tri.tmp" diff-${testname}.out diff-${testname}.out diff-${label}+ext_hfirst-${iext_hfirst} ""
      else
         petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
      fi

   done
done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
