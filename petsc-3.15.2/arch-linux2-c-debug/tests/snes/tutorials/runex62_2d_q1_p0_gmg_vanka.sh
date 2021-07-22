#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex62'
testname='runex62_2d_q1_p0_gmg_vanka'
label='snes_tutorials-ex62_2d_q1_p0_gmg_vanka'
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
args='-sol quadratic -dm_plex_box_simplex 0 -dm_refine_hierarchy 2 -vel_petscspace_degree 1 -pres_petscspace_degree 0 -petscds_jac_pre 0 -snes_rtol 1.0e-4 -ksp_type fgmres -ksp_atol 1e-5 -ksp_error_if_not_converged -pc_type mg -mg_levels_ksp_type gmres -mg_levels_ksp_max_it 30 -mg_levels_pc_type patch -mg_levels_pc_patch_partition_of_unity 0 -mg_levels_pc_patch_construct_codim 0 -mg_levels_pc_patch_construct_type vanka -mg_levels_sub_ksp_type preonly -mg_levels_sub_pc_type lu -mg_coarse_pc_type svd'
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

   petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args} " ex62_2d_q1_p0_gmg_vanka.tmp ${testname}.err "${label}" 
   res=$?

   if test $res = 0; then
      petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/snes/tutorials/output/ex62_2d_q1_p0_gmg_vanka.out ex62_2d_q1_p0_gmg_vanka.tmp" diff-${testname}.out diff-${testname}.out diff-${label} ""
   else
      petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
   fi

done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
