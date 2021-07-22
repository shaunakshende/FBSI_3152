#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex13'
testname='runex13_aijmkl_seq'
label='snes_tests-ex13_aijmkl_seq'
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
args='-dm_plex_box_dim 3 -dm_plex_box_faces 4,4,4 -dm_refine 1 -petscpartitioner_type simple -potential_petscspace_degree 1 -dm_distribute -dm_plex_box_simplex 0 -snes_monitor_short -snes_type ksponly -dm_view -pc_type gamg -pc_gamg_sym_graph 0 -pc_gamg_threshold -1 -pc_gamg_square_graph 10 -pc_gamg_process_eq_limit 400 -pc_gamg_reuse_interpolation -pc_gamg_coarse_eq_limit 10 -mg_levels_esteig_ksp_type cg -mg_levels_pc_type jacobi -ksp_type cg -ksp_norm_type unpreconditioned -snes_converged_reason -ksp_converged_reason -snes_rtol 1.e-4 -dm_mat_type aijmkl -dm_vec_type standard'
diff_args=''
timeoutfactor=1

mpiexec=${PETSCMPIEXEC:-"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/bin/mpiexec"}
diffexec=${PETSCDIFF:-"${petsc_bindir}/petscdiff"}

. "${config_dir}/petsc_harness.sh"

# The diff flags come from script arguments
diff_exe="${diffexec} ${diff_flags} ${diff_args}"
mpiexec="${mpiexec} ${mpiexec_flags}"

if ! $force; then
    petsc_report_tapoutput "" "${label}" "TODO broken (INDEFINITE PC)"
    total=1; todo=1
    petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
    exit
fi


if ! $force; then
    petsc_report_tapoutput "" "${label}" "SKIP PETSC_HAVE_MKL_SPARSE requirement not met"
    total=1; skip=1
    petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
    exit
fi


nsize_in=${nsize:-1}


for insize in ${nsize_in}; do

   petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args} " ex13_aijmkl_seq.tmp ${testname}.err "${label}" 
   res=$?

   if test $res = 0; then
      petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/snes/tests/output/ex13_aijmkl_seq.out ex13_aijmkl_seq.tmp" diff-${testname}.out diff-${testname}.out diff-${label} ""
   else
      petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
   fi

done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
