#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex56'
testname='runex56_viennacl'
label='snes_tutorials-ex56_viennacl'
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
args='-cells 2,2,1 -max_conv_its 2 -petscspace_degree 2 -snes_max_it 2 -ksp_max_it 100 -ksp_type cg -ksp_rtol 1.e-10 -ksp_norm_type unpreconditioned -snes_rtol 1.e-10 -pc_type gamg -pc_gamg_type agg -pc_gamg_agg_nsmooths 1 -pc_gamg_coarse_eq_limit 10 -pc_gamg_reuse_interpolation true -pc_gamg_square_graph 1 -pc_gamg_threshold 0.05 -pc_gamg_threshold_scale .0 -use_mat_nearnullspace true -mg_levels_ksp_max_it 2 -mg_levels_ksp_type chebyshev -mg_levels_ksp_chebyshev_esteig 0,0.05,0,1.05 -mg_levels_pc_type jacobi -ex56_dm_mat_type aijviennacl -ex56_dm_vec_type viennacl -ksp_monitor_short -ksp_converged_reason -snes_converged_reason -snes_monitor_short -ex56_dm_view -petscpartitioner_type simple -pc_gamg_process_eq_limit 20'
diff_args=''
timeoutfactor=1

mpiexec=${PETSCMPIEXEC:-"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/bin/mpiexec"}
diffexec=${PETSCDIFF:-"${petsc_bindir}/petscdiff"}

. "${config_dir}/petsc_harness.sh"

# The diff flags come from script arguments
diff_exe="${diffexec} ${diff_flags} ${diff_args}"
mpiexec="${mpiexec} ${mpiexec_flags}"



if ! $force; then
    petsc_report_tapoutput "" "${label}" "SKIP PETSC_HAVE_VIENNACL requirement not met"
    total=1; skip=1
    petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
    exit
fi


nsize_in=${nsize:-4}


for insize in ${nsize_in}; do

   petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args} " ex56_viennacl.tmp ${testname}.err "${label}" 
   res=$?

   if test $res = 0; then
      petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/snes/tutorials/output/ex56_viennacl.out ex56_viennacl.tmp" diff-${testname}.out diff-${testname}.out diff-${label} ""
   else
      petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
   fi

done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
