#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex69'
testname='runex69_p2p1fetidp_olof_additive'
label='snes_tutorials-ex69_p2p1fetidp_olof_additive'
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
args='-dm_distribute -dm_plex_separate_marker -dm_refine_pre 2 -dm_refine_uniform_pre -dm_mat_type is -dm_view -petscpartitioner_type simple -vel_petscspace_degree 2 -pres_petscspace_degree 1 -snes_error_if_not_converged -ksp_error_if_not_converged -ksp_type fetidp -ksp_fetidp_saddlepoint -fetidp_ksp_type cg -fetidp_ksp_norm_type natural -fetidp_pc_discrete_harmonic 1 -fetidp_harmonic_pc_type cholesky -ksp_fetidp_pressure_schur -fetidp_fieldsplit_p_pc_type bddc -fetidp_fieldsplit_p_pc_bddc_dirichlet_pc_type none -fetidp_fieldsplit_p_ksp_type preonly -fetidp_fieldsplit_lag_ksp_error_if_not_converged 0 -fetidp_fieldsplit_lag_ksp_type chebyshev -fetidp_fieldsplit_lag_ksp_max_it 2 -fetidp_bddc_pc_bddc_detect_disconnected -fetidp_bddc_pc_bddc_symmetric -fetidp_bddc_pc_bddc_vertex_size 3 -fetidp_bddc_pc_bddc_graph_maxcount 2 -fetidp_bddc_pc_bddc_coarse_redundant_pc_type cholesky -fetidp_bddc_pc_bddc_dirichlet_pc_type none -fetidp_bddc_pc_bddc_neumann_pc_type svd -fetidp_pc_fieldsplit_type additive'
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


nsize_in=${nsize:-4}


for insize in ${nsize_in}; do

   petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args} " ex69_p2p1fetidp_olof_additive.tmp ${testname}.err "${label}" 
   res=$?

   if test $res = 0; then
      petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/snes/tutorials/output/ex69_p2p1fetidp_olof.out ex69_p2p1fetidp_olof_additive.tmp" diff-${testname}.out diff-${testname}.out diff-${label} ""
   else
      petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
   fi

done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
