#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex1f90'
testname='runex1f90_0'
label='ts_utils_dmplexlandau_tutorials-ex1f90_0'
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
args='-petscspace_degree 3 -petscspace_poly_tensor 1 -dm_landau_type p4est -info :dm,tsadapt -dm_landau_ion_masses 2,4 -dm_landau_ion_charges 1,18 -dm_landau_thermal_temps 5,5,.5 -dm_landau_n 1.00018,1,1e-5 -dm_landau_n_0 1e20 -ex1f90_ts_monitor -ex1f90_snes_rtol 1.e-6 -ex1f90_snes_monitor -ex1f90_snes_converged_reason -ex1f90_ts_type arkimex -ex1f90_ts_arkimex_type 1bee -ex1f90_ts_max_snes_failures -1 -ex1f90_ts_rtol 1e-4 -ex1f90_ts_dt 1.e-6 -ex1f90_ts_max_time 1 -ex1f90_ts_adapt_clip .5,1.25 -ex1f90_ts_adapt_scale_solve_failed 0.75 -ex1f90_ts_adapt_time_step_increase_delay 5 -ex1f90_ts_max_steps 1 -ex1f90_pc_type lu -ex1f90_ksp_type preonly -dm_landau_amr_levels_max 7 -dm_landau_domain_radius 5 -dm_landau_amr_re_levels 0 -dm_landau_re_radius 1 -dm_landau_amr_z_refine1 1 -dm_landau_amr_z_refine2 0 -dm_landau_amr_post_refine 0 -dm_landau_z_radius1 .1 -dm_landau_z_radius2 .1 -dm_refine 1 -dm_landau_device_type cpu'
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


nsize_in=${nsize:-1}


for insize in ${nsize_in}; do

   petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args} " ex1f90_0.tmp ${testname}.err "${label}" 
   res=$?

   if test $res = 0; then
      petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/ts/utils/dmplexlandau/tutorials/output/ex1f90_0.out ex1f90_0.tmp" diff-${testname}.out diff-${testname}.out diff-${label} ""
   else
      petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
   fi

done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
