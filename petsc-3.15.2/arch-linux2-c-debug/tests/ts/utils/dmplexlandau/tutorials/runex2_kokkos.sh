#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex2'
testname='runex2_kokkos'
label='ts_utils_dmplexlandau_tutorials-ex2_kokkos'
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
args='-dm_landau_Ez 0 -petscspace_degree 2 -petscspace_poly_tensor 1 -dm_landau_type p4est -info :dm,tsadapt -dm_landau_ion_masses 2 -dm_landau_ion_charges 1 -dm_landau_thermal_temps 5,5 -dm_landau_n 2,2 -dm_landau_n_0 5e19 -ts_monitor -snes_rtol 1.e-10 -snes_stol 1.e-14 -snes_monitor -snes_converged_reason -snes_max_it 10 -ts_type arkimex -ts_arkimex_type 1bee -ts_max_snes_failures -1 -ts_rtol 1e-3 -ts_dt 1.e-1 -ts_max_time 1 -ts_adapt_clip .5,1.25 -ts_max_steps 2 -ts_adapt_scale_solve_failed 0.75 -ts_adapt_time_step_increase_delay 5 -pc_type lu -ksp_type preonly -dm_landau_amr_levels_max 9 -dm_landau_domain_radius -.75 -ex2_impurity_source_type pulse -ex2_pulse_start_time 1e-1 -ex2_pulse_width_time 10 -ex2_pulse_rate 1e-2 -ex2_t_cold .05 -ex2_plot_dt 1e-1 -dm_refine 1 -dm_preallocate_only -dm_landau_device_type kokkos -dm_landau_gpu_assembly true -dm_mat_type aijkokkos -dm_vec_type kokkos'
diff_args=''
timeoutfactor=1

mpiexec=${PETSCMPIEXEC:-"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/bin/mpiexec"}
diffexec=${PETSCDIFF:-"${petsc_bindir}/petscdiff"}

. "${config_dir}/petsc_harness.sh"

# The diff flags come from script arguments
diff_exe="${diffexec} ${diff_flags} ${diff_args}"
mpiexec="${mpiexec} ${mpiexec_flags}"



if ! $force; then
    petsc_report_tapoutput "" "${label}" "SKIP PETSC_HAVE_P4EST requirement not met, PETSC_HAVE_KOKKOS_KERNELS requirement not met, Null requirement not met: define(PETSC_USE_CTABLE)"
    total=1; skip=1
    petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
    exit
fi


nsize_in=${nsize:-1}


for insize in ${nsize_in}; do

   petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args} " ex2_kokkos.tmp ${testname}.err "${label}" 
   res=$?

   if test $res = 0; then
      petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/ts/utils/dmplexlandau/tutorials/output/ex2_kokkos.out ex2_kokkos.tmp" diff-${testname}.out diff-${testname}.out diff-${label} ""
   else
      petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
   fi

done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
