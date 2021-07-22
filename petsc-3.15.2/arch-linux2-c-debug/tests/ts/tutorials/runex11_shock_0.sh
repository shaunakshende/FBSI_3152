#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex11'
testname='runex11_shock_0'
label='ts_tutorials-ex11_shock_0'
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
args='-ufv_vtk_interval 0 -monitor density,energy -f -grid_size 2,1 -grid_bounds -1,1.,0.,1 -bc_wall 1,2,3,4 -dm_type p4est -dm_forest_partition_overlap 1 -dm_forest_maximum_refinement 6 -dm_forest_minimum_refinement 2 -dm_forest_initial_refinement 2 -ufv_use_amr -refine_vec_tagger_box 0.5,inf -coarsen_vec_tagger_box 0,1.e-2 -refine_tag_view -coarsen_tag_view -physics euler -eu_type iv_shock -ufv_cfl 10 -eu_alpha 60. -grid_skew_60 -eu_gamma 1.4 -eu_amach 2.02 -eu_rho2 3. -petscfv_type leastsquares -petsclimiter_type minmod -petscfv_compute_gradients 0 -ts_max_time 0.5 -ts_ssp_type rks2 -ts_ssp_nstages 10 -ufv_vtk_basename ${wPETSC_DIR}/ex11'
diff_args=''
timeoutfactor=3

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


nsize_in=${nsize:-1}


for insize in ${nsize_in}; do

   petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args} " ex11_shock_0.tmp ${testname}.err "${label}" 
   res=$?

   if test $res = 0; then
      petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/ts/tutorials/output/ex11_shock_0.out ex11_shock_0.tmp" diff-${testname}.out diff-${testname}.out diff-${label} ""
   else
      petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
   fi

done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
