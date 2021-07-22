#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex62'
testname='runex62_3d_p2_p1_fetidp'
label='snes_tutorials-ex62_3d_p2_p1_fetidp'
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
args='-sol quadratic -dm_plex_box_dim 3 -dm_plex_box_faces 2,2,2 -dm_refine 1 -dm_mat_type is -dm_distribute -petscpartitioner_type simple -vel_petscspace_degree 2 -pres_petscspace_degree 1 -petscds_jac_pre 0 -snes_error_if_not_converged -ksp_type fetidp -ksp_rtol 1.0e-9 -ksp_error_if_not_converged -ksp_fetidp_saddlepoint -fetidp_ksp_type cg -fetidp_fieldsplit_p_ksp_max_it 1 -fetidp_fieldsplit_p_ksp_type richardson -fetidp_fieldsplit_p_ksp_richardson_scale 1000 -fetidp_fieldsplit_p_pc_type none -fetidp_bddc_pc_bddc_use_deluxe_scaling -fetidp_bddc_pc_bddc_benign_trick -fetidp_bddc_pc_bddc_deluxe_singlemat -fetidp_pc_discrete_harmonic -fetidp_harmonic_pc_factor_mat_solver_type petsc -fetidp_harmonic_pc_type cholesky -fetidp_bddelta_pc_factor_mat_solver_type umfpack -fetidp_fieldsplit_lag_ksp_type preonly -fetidp_bddc_sub_schurs_mat_solver_type mumps -fetidp_bddc_sub_schurs_mat_mumps_icntl_14 100000 -fetidp_bddelta_pc_factor_mat_ordering_type external -fetidp_bddc_pc_bddc_dirichlet_pc_factor_mat_solver_type umfpack -fetidp_bddc_pc_bddc_neumann_pc_factor_mat_solver_type umfpack -fetidp_bddc_pc_bddc_dirichlet_pc_factor_mat_ordering_type external -fetidp_bddc_pc_bddc_neumann_pc_factor_mat_ordering_type external'
diff_args=''
timeoutfactor=1

mpiexec=${PETSCMPIEXEC:-"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/bin/mpiexec"}
diffexec=${PETSCDIFF:-"${petsc_bindir}/petscdiff"}

. "${config_dir}/petsc_harness.sh"

# The diff flags come from script arguments
diff_exe="${diffexec} ${diff_flags} ${diff_args}"
mpiexec="${mpiexec} ${mpiexec_flags}"



if ! $force; then
    petsc_report_tapoutput "" "${label}" "SKIP PETSC_HAVE_CTETGEN requirement not met, PETSC_HAVE_MUMPS requirement not met, PETSC_HAVE_SUITESPARSE requirement not met"
    total=1; skip=1
    petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
    exit
fi


nsize_in=${nsize:-5}


for insize in ${nsize_in}; do

   petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args} " ex62_3d_p2_p1_fetidp.tmp ${testname}.err "${label}" 
   res=$?

   if test $res = 0; then
      petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/snes/tutorials/output/ex62_3d_p2_p1_fetidp.out ex62_3d_p2_p1_fetidp.tmp" diff-${testname}.out diff-${testname}.out diff-${label} ""
   else
      petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
   fi

done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
