#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex12'
testname='runex12_p4est_exact_q2_nonconformal_parallel'
label='snes_tutorials-ex12_p4est_exact_q2_nonconformal_parallel'
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
args='-run_type exact -interpolate 1 -petscspace_degree 2 -fas_levels_snes_atol 1.e-10 -snes_max_it 1 -snes_type fas -snes_fas_levels 3 -fas_coarse_pc_type none -fas_coarse_ksp_type preonly -fas_coarse_snes_monitor_short -snes_monitor_short -fas_coarse_snes_linesearch_type basic -snes_converged_reason ::ascii_info_detail -snes_view -fas_levels_snes_type newtonls -fas_levels_pc_type none -fas_levels_ksp_type preonly -fas_levels_snes_monitor_short -simplex 0 -dm_plex_convert_type p4est -dm_forest_minimum_refinement 0 -dm_forest_initial_refinement 2 -dm_forest_maximum_refinement 4 -dm_p4est_refine_pattern hash -petscpartitioner_type simple -cells 2,2'
diff_args=''
timeoutfactor=1

mpiexec=${PETSCMPIEXEC:-"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/bin/mpiexec"}
diffexec=${PETSCDIFF:-"${petsc_bindir}/petscdiff"}

. "${config_dir}/petsc_harness.sh"

# The diff flags come from script arguments
diff_exe="${diffexec} ${diff_flags} ${diff_args}"
mpiexec="${mpiexec} ${mpiexec_flags}"



if ! $force; then
    petsc_report_tapoutput "" "${label}" "SKIP PETSC_HAVE_P4EST requirement not met"
    total=1; skip=1
    petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
    exit
fi


nsize_in=${nsize:-7}


for insize in ${nsize_in}; do

   petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args} " ex12_p4est_exact_q2_nonconformal_parallel.tmp ${testname}.err "${label}" 
   res=$?

   if test $res = 0; then
      petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/snes/tutorials/output/ex12_p4est_exact_q2_nonconformal_parallel.out ex12_p4est_exact_q2_nonconformal_parallel.tmp" diff-${testname}.out diff-${testname}.out diff-${label} ""
   else
      petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
   fi

done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
