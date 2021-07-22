#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex62'
testname='runex62_3d_q1_p0_conv'
label='snes_tutorials-ex62_3d_q1_p0_conv'
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
args='-sol quadratic -dm_plex_box_simplex 0 -dm_plex_box_dim 3 -dm_refine 1 -vel_petscspace_degree 1 -pres_petscspace_degree 0 -snes_convergence_estimate -convest_num_refine 1 -ksp_atol 1e-10 -ksp_error_if_not_converged -petscds_jac_pre 0 -pc_type fieldsplit -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type full -pc_fieldsplit_schur_precondition full -fieldsplit_velocity_pc_type lu -fieldsplit_pressure_ksp_rtol 1e-10 -fieldsplit_pressure_pc_type gamg -fieldsplit_pressure_mg_levels_pc_type jacobi -fieldsplit_pressure_mg_coarse_pc_type svd'
diff_args=''
timeoutfactor=1

mpiexec=${PETSCMPIEXEC:-"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/bin/mpiexec"}
diffexec=${PETSCDIFF:-"${petsc_bindir}/petscdiff"}

. "${config_dir}/petsc_harness.sh"

# The diff flags come from script arguments
diff_exe="${diffexec} ${diff_flags} ${diff_args}"
mpiexec="${mpiexec} ${mpiexec_flags}"





nsize_in=${nsize:-1}


for insize in ${nsize_in}; do

   petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args} " ex62_3d_q1_p0_conv.tmp ${testname}.err "${label}" 
   res=$?

   if test $res = 0; then
      petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/snes/tutorials/output/ex62_3d_q1_p0_conv.out ex62_3d_q1_p0_conv.tmp" diff-${testname}.out diff-${testname}.out diff-${label} ""
   else
      petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
   fi

done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
