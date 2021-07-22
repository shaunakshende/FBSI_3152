#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex39'
testname='runex39_quad_rt_0'
label='dm_impls_plex_tests-ex39_quad_rt_0'
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
args='-dim 2 -simplex false -mesh_transform skew -divErr_petscspace_degree 1 -divErr_petscdualspace_lagrange_continuity false -dm_refine 0 -snes_error_if_not_converged -ksp_rtol 1e-10 -ksp_error_if_not_converged -pc_type fieldsplit -pc_fieldsplit_detect_saddle_point -pc_fieldsplit_type schur -pc_fieldsplit_schur_precondition full -velocity_petscfe_default_quadrature_order 1 -velocity_petscspace_type sum -velocity_petscspace_variables 2 -velocity_petscspace_components 2 -velocity_petscspace_sum_spaces 2 -velocity_petscspace_sum_concatenate true -velocity_subspace0_petscspace_variables 2 -velocity_subspace0_petscspace_type tensor -velocity_subspace0_petscspace_tensor_spaces 2 -velocity_subspace0_petscspace_tensor_uniform false -velocity_subspace0_subspace_0_petscspace_degree 1 -velocity_subspace0_subspace_1_petscspace_degree 0 -velocity_subspace1_petscspace_variables 2 -velocity_subspace1_petscspace_type tensor -velocity_subspace1_petscspace_tensor_spaces 2 -velocity_subspace1_petscspace_tensor_uniform false -velocity_subspace1_subspace_0_petscspace_degree 0 -velocity_subspace1_subspace_1_petscspace_degree 1 -velocity_petscdualspace_form_degree -1 -velocity_petscdualspace_order 1 -velocity_petscdualspace_lagrange_trimmed true'
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

   petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args} " ex39_quad_rt_0.tmp ${testname}.err "${label}" 
   res=$?

   if test $res = 0; then
      petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/dm/impls/plex/tests/output/ex39_quad_rt_0.out ex39_quad_rt_0.tmp" diff-${testname}.out diff-${testname}.out diff-${label} ""
   else
      petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
   fi

done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
