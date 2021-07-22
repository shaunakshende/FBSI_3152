#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex76'
testname='runex76_geneo_share_cholesky'
label='ksp_ksp_tutorials-ex76_geneo_share_cholesky'
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
args='-ksp_converged_reason -ksp_max_it 150 -pc_type hpddm -pc_hpddm_levels_1_eps_nev 5 -pc_hpddm_coarse_p 1 -pc_hpddm_coarse_pc_type redundant -load_dir ${DATAFILESPATH}/matrices/hpddm/GENEO -pc_hpddm_define_subdomains -pc_hpddm_has_neumann -pc_hpddm_levels_1_sub_pc_type cholesky -pc_hpddm_levels_1_st_pc_type cholesky -pc_hpddm_levels_1_eps_gen_non_hermitian'
diff_args=''
timeoutfactor=1

mpiexec=${PETSCMPIEXEC:-"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/bin/mpiexec"}
diffexec=${PETSCDIFF:-"${petsc_bindir}/petscdiff"}

. "${config_dir}/petsc_harness.sh"

# The diff flags come from script arguments
diff_exe="${diffexec} ${diff_flags} ${diff_args}"
mpiexec="${mpiexec} ${mpiexec_flags}"



if ! $force; then
    petsc_report_tapoutput "" "${label}" "SKIP PETSC_HAVE_HPDDM requirement not met, PETSC_HAVE_SLEPC requirement not met, Requires DATAFILESPATH"
    total=1; skip=1
    petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
    exit
fi


nsize_in=${nsize:-4}
mat_type_in=${mat_type:-aij baij sbaij}
pc_hpddm_levels_1_st_share_sub_ksp_in=${pc_hpddm_levels_1_st_share_sub_ksp:-false true}


for insize in ${nsize_in}; do
   for imat_type in ${mat_type_in}; do
      for ipc_hpddm_levels_1_st_share_sub_ksp in ${pc_hpddm_levels_1_st_share_sub_ksp_in}; do

         petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args}  -mat_type ${imat_type} -pc_hpddm_levels_1_st_share_sub_ksp ${ipc_hpddm_levels_1_st_share_sub_ksp}" ex76_geneo_share_cholesky.tmp ${testname}.err "${label}+mat_type-${imat_type}_pc_hpddm_levels_1_st_share_sub_ksp-${ipc_hpddm_levels_1_st_share_sub_ksp}" 
         res=$?

         if test $res = 0; then
            petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/ksp/ksp/tutorials/output/ex76_geneo_share.out ex76_geneo_share_cholesky.tmp" diff-${testname}.out diff-${testname}.out diff-${label}+mat_type-${imat_type}_pc_hpddm_levels_1_st_share_sub_ksp-${ipc_hpddm_levels_1_st_share_sub_ksp} ""
         else
            petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
         fi

      done
   done
done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
