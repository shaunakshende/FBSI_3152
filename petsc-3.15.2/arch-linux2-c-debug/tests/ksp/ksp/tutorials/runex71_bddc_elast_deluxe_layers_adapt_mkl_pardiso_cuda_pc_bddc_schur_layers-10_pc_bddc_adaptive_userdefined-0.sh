#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex71'
testname='runex71_bddc_elast_deluxe_layers_adapt_mkl_pardiso_cuda_pc_bddc_schur_layers-10_pc_bddc_adaptive_userdefined-0'
label='ksp_ksp_tutorials-ex71_bddc_elast_deluxe_layers_adapt_mkl_pardiso_cuda_pc_bddc_schur_layers-10_pc_bddc_adaptive_userdefined-0'
runfiles=''
wPETSC_DIR='/home/shaunak/Desktop/petsc-3.15.2'
petsc_dir='/home/shaunak/Desktop/petsc-3.15.2'
petsc_arch='arch-linux2-c-debug'
# Must be consistent with gmakefile.test
testlogtapfile=/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests/test_${petsc_arch}_tap.log
testlogerrfile=/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests/test_${petsc_arch}_err.log
config_dir='/home/shaunak/Desktop/petsc-3.15.2/config'
filter='sed -e "s/CONVERGED_RTOL iterations 6/CONVERGED_RTOL iterations 5/g"'
filter_output=''
petsc_bindir='/home/shaunak/Desktop/petsc-3.15.2/lib/petsc/bin'
DATAFILESPATH=${DATAFILESPATH:-""}
args='-pde_type Elasticity -cells 7,9,8 -dim 3 -ksp_converged_reason -pc_bddc_coarse_redundant_pc_type svd -ksp_error_if_not_converged -pc_bddc_monolithic -sub_schurs_mat_solver_type mkl_pardiso -sub_schurs_mat_mkl_pardiso_65 1 -pc_bddc_use_deluxe_scaling -pc_bddc_adaptive_threshold 2.0 -matis_localmat_type seqaijcusparse -pc_bddc_schur_layers 10 -pc_bddc_adaptive_userdefined 0'
diff_args=''
timeoutfactor=1

mpiexec=${PETSCMPIEXEC:-"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/bin/mpiexec"}
diffexec=${PETSCDIFF:-"${petsc_bindir}/petscdiff"}

. "${config_dir}/petsc_harness.sh"

# The diff flags come from script arguments
diff_exe="${diffexec} ${diff_flags} ${diff_args}"
mpiexec="${mpiexec} ${mpiexec_flags}"



if ! $force; then
    petsc_report_tapoutput "" "${label}" "SKIP PETSC_HAVE_MKL_PARDISO requirement not met, PETSC_HAVE_CUDA requirement not met, Required: define(PETSC_HAVE_CUSOLVERDNDPOTRI)"
    total=1; skip=1
    petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
    exit
fi


nsize_in=${nsize:-8}
sub_schurs_schur_mat_type_in=${sub_schurs_schur_mat_type:-seqdensecuda seqdense}


for insize in ${nsize_in}; do
   for isub_schurs_schur_mat_type in ${sub_schurs_schur_mat_type_in}; do

      petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args}  -sub_schurs_schur_mat_type ${isub_schurs_schur_mat_type}" ex71_bddc_elast_deluxe_layers_adapt_mkl_pardiso_cuda_pc_bddc_schur_layers-10_pc_bddc_adaptive_userdefined-0.tmp ${testname}.err "${label}+sub_schurs_schur_mat_type-${isub_schurs_schur_mat_type}" 
      res=$?

      if test $res = 0; then
         petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/ksp/ksp/tutorials/output/ex71_bddc_elast_deluxe_layers_adapt_mkl_pardiso_cuda_pc_bddc_schur_layers-10_pc_bddc_adaptive_userdefined-0.out ex71_bddc_elast_deluxe_layers_adapt_mkl_pardiso_cuda_pc_bddc_schur_layers-10_pc_bddc_adaptive_userdefined-0.tmp" diff-${testname}.out diff-${testname}.out diff-${label}+sub_schurs_schur_mat_type-${isub_schurs_schur_mat_type} ""
      else
         petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
      fi

   done
done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
