#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex59'
testname='runex59_bddc_fetidp_ml_eqlimit_1'
label='ksp_ksp_tutorials-ex59_bddc_fetidp_ml_eqlimit_1'
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
args='-npx 3 -npy 3 -p 2 -nex 6 -ney 6 -physical_pc_bddc_deluxe_singlemat -physical_sub_schurs_mat_solver_type petsc -physical_pc_bddc_use_deluxe_scaling -physical_pc_bddc_graph_maxcount 1 -physical_pc_bddc_levels 3 -physical_pc_bddc_coarsening_ratio 2 -physical_pc_bddc_coarse_ksp_type gmres -physical_pc_bddc_coarse_eqs_limit 31 -mat_partitioning_type average -physical_pc_bddc_coarse_pc_bddc_graph_maxcount 1'
diff_args=''
timeoutfactor=1

mpiexec=${PETSCMPIEXEC:-"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/bin/mpiexec"}
diffexec=${PETSCDIFF:-"${petsc_bindir}/petscdiff"}

. "${config_dir}/petsc_harness.sh"

# The diff flags come from script arguments
diff_exe="${diffexec} ${diff_flags} ${diff_args}"
mpiexec="${mpiexec} ${mpiexec_flags}"





nsize_in=${nsize:-9}


for insize in ${nsize_in}; do

   petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args} " ex59_bddc_fetidp_ml_eqlimit_1.tmp ${testname}.err "${label}" 
   res=$?

   if test $res = 0; then
      petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/ksp/ksp/tutorials/output/ex59_bddc_fetidp_ml_eqlimit.out ex59_bddc_fetidp_ml_eqlimit_1.tmp" diff-${testname}.out diff-${testname}.out diff-${label} ""
   else
      petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
   fi

done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 