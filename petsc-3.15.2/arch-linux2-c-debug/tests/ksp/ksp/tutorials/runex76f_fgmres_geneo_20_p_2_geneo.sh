#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex76f'
testname='runex76f_fgmres_geneo_20_p_2_geneo'
label='ksp_ksp_tutorials-ex76f_fgmres_geneo_20_p_2_geneo'
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
args='-ksp_converged_reason -pc_type hpddm -pc_hpddm_levels_1_sub_pc_type cholesky -pc_hpddm_levels_1_eps_nev 20 -pc_hpddm_levels_2_p 2 -pc_hpddm_levels_2_sub_pc_type cholesky -pc_hpddm_levels_2_ksp_type gmres -ksp_type fgmres -load_dir ${DATAFILESPATH}/matrices/hpddm/GENEO'
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
pc_hpddm_levels_2_mat_type_in=${pc_hpddm_levels_2_mat_type:-baij sbaij}
pc_hpddm_levels_2_eps_nev_in=${pc_hpddm_levels_2_eps_nev:-5 20}
pc_hpddm_coarse_mat_type_in=${pc_hpddm_coarse_mat_type:-baij sbaij}


for insize in ${nsize_in}; do
   for ipc_hpddm_levels_2_mat_type in ${pc_hpddm_levels_2_mat_type_in}; do
      for ipc_hpddm_levels_2_eps_nev in ${pc_hpddm_levels_2_eps_nev_in}; do
         for ipc_hpddm_coarse_mat_type in ${pc_hpddm_coarse_mat_type_in}; do

            petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args}  -pc_hpddm_levels_2_mat_type ${ipc_hpddm_levels_2_mat_type} -pc_hpddm_levels_2_eps_nev ${ipc_hpddm_levels_2_eps_nev} -pc_hpddm_coarse_mat_type ${ipc_hpddm_coarse_mat_type}" ex76f_fgmres_geneo_20_p_2_geneo.tmp ${testname}.err "${label}+pc_hpddm_levels_2_mat_type-${ipc_hpddm_levels_2_mat_type}_pc_hpddm_levels_2_eps_nev-${ipc_hpddm_levels_2_eps_nev}_pc_hpddm_coarse_mat_type-${ipc_hpddm_coarse_mat_type}" 
            res=$?

            if test $res = 0; then
               petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/ksp/ksp/tutorials/output/ex76_fgmres_geneo_20_p_2.out ex76f_fgmres_geneo_20_p_2_geneo.tmp" diff-${testname}.out diff-${testname}.out diff-${label}+pc_hpddm_levels_2_mat_type-${ipc_hpddm_levels_2_mat_type}_pc_hpddm_levels_2_eps_nev-${ipc_hpddm_levels_2_eps_nev}_pc_hpddm_coarse_mat_type-${ipc_hpddm_coarse_mat_type} ""
            else
               petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
            fi

         done
      done
   done
done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
