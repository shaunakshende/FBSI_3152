#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex69'
testname='runex69_q2p1bddc_benign_deluxe_mkl'
label='snes_tutorials-ex69_q2p1bddc_benign_deluxe_mkl'
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
args='-dm_distribute -dm_plex_box_simplex 0 -dm_plex_separate_marker -dm_refine_pre 1 -dm_mat_type is -dm_view -petscpartitioner_type simple -vel_petscspace_degree 2 -pres_petscspace_degree 1 -pres_petscspace_poly_tensor 0 -pres_petscdualspace_lagrange_continuity 0 -pres_petscdualspace_lagrange_node_endpoints 0 -petscds_jac_pre 0 -snes_error_if_not_converged -ksp_error_if_not_converged -ksp_type cg -ksp_norm_type natural -pc_type bddc -pc_bddc_benign_trick -pc_bddc_nonetflux -pc_bddc_detect_disconnected -pc_bddc_vertex_size 2 -pc_bddc_coarse_redundant_pc_type svd -pc_bddc_use_qr_single -pc_bddc_use_deluxe_scaling -pc_bddc_deluxe_zerorows -sub_schurs_mat_solver_type mkl_pardiso -snes_rtol 1.e-7'
diff_args=''
timeoutfactor=1

mpiexec=${PETSCMPIEXEC:-"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/bin/mpiexec"}
diffexec=${PETSCDIFF:-"${petsc_bindir}/petscdiff"}

. "${config_dir}/petsc_harness.sh"

# The diff flags come from script arguments
diff_exe="${diffexec} ${diff_flags} ${diff_args}"
mpiexec="${mpiexec} ${mpiexec_flags}"



if ! $force; then
    petsc_report_tapoutput "" "${label}" "SKIP PETSC_HAVE_MKL_PARDISO requirement not met"
    total=1; skip=1
    petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
    exit
fi


nsize_in=${nsize:-5}


for insize in ${nsize_in}; do

   petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args} " ex69_q2p1bddc_benign_deluxe_mkl.tmp ${testname}.err "${label}" 
   res=$?

   if test $res = 0; then
      petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/snes/tutorials/output/ex69_q2p1fetidp_deluxe.out ex69_q2p1bddc_benign_deluxe_mkl.tmp" diff-${testname}.out diff-${testname}.out diff-${label} ""
   else
      petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
   fi

done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
