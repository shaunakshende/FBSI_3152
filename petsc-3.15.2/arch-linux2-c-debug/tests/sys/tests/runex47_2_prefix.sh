#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex47'
testname='runex47_2_prefix'
label='sys_tests-ex47_2_prefix'
runfiles='ex47-opt.txt ex47-opt.yml ex47-opt.json'
wPETSC_DIR='/home/shaunak/Desktop/petsc-3.15.2'
petsc_dir='/home/shaunak/Desktop/petsc-3.15.2'
petsc_arch='arch-linux2-c-debug'
# Must be consistent with gmakefile.test
testlogtapfile=/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests/test_${petsc_arch}_tap.log
testlogerrfile=/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests/test_${petsc_arch}_err.log
config_dir='/home/shaunak/Desktop/petsc-3.15.2/config'
filter='egrep -v "(options_left|options_monitor|malloc_dump|malloc_test|saws_port_auto_select|display|check_pointer_intensity|error_output_stdout|nox|vecscatter_mpi1|use_gpu_aware_mpi)"'
filter_output=''
petsc_bindir='/home/shaunak/Desktop/petsc-3.15.2/lib/petsc/bin'
DATAFILESPATH=${DATAFILESPATH:-""}
args='-options_left false -options_monitor -options_file ex47-opt.txt -prefix_push p5_ -options_file ex47-opt.yml -prefix_pop -prefix_push p5_ -options_file ex47-opt.yml:yaml -prefix_pop -prefix_push p6_ -options_file_yaml ex47-opt.yml -prefix_pop -prefix_push p7_ -options_string_yaml "`cat ex47-opt.yml`" -prefix_pop -prefix_push p7_ -options_string_yaml "`cat ex47-opt.yml`" -prefix_pop -prefix_push p8_ -options_string_yaml "`cat ex47-opt.yml`" -prefix_pop -prefix_push p9_ -options_file ex47-opt.json -prefix_pop'
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

   petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args} " ex47_2_prefix.tmp ${testname}.err "${label}" 
   res=$?

   if test $res = 0; then
      petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/sys/tests/output/ex47_2_prefix.out ex47_2_prefix.tmp" diff-${testname}.out diff-${testname}.out diff-${label} ""
   else
      petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
   fi

done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 