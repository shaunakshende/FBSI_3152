#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex76'
testname='runex76_2d_tri_p2_p1_p1_conv'
label='snes_tutorials-ex76_2d_tri_p2_p1_p1_conv'
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
args='-dm_plex_separate_marker -sol_type cubic -dm_refine 0 -vel_petscspace_degree 2 -pres_petscspace_degree 1 -temp_petscspace_degree 1 -snes_error_if_not_converged -snes_convergence_test correct_pressure -snes_convergence_estimate -convest_num_refine 1 -ksp_type fgmres -ksp_gmres_restart 10 -ksp_rtol 1.0e-9 -ksp_error_if_not_converged -pc_type fieldsplit -pc_fieldsplit_0_fields 0,2 -pc_fieldsplit_1_fields 1 -pc_fieldsplit_type schur -pc_fieldsplit_schur_factorization_type full -fieldsplit_0_pc_type lu -fieldsplit_pressure_ksp_rtol 1e-10 -fieldsplit_pressure_pc_type jacobi'
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


nsize_in=${nsize:-1}


for insize in ${nsize_in}; do

   petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args} " ex76_2d_tri_p2_p1_p1_conv.tmp ${testname}.err "${label}" 
   res=$?

   if test $res = 0; then
      petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/snes/tutorials/output/ex76_2d_tri_p2_p1_p1_conv.out ex76_2d_tri_p2_p1_p1_conv.tmp" diff-${testname}.out diff-${testname}.out diff-${label} ""
   else
      petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
   fi

done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
