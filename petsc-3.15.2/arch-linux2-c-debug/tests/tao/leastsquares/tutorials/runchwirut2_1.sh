#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../chwirut2'
testname='runchwirut2_1'
label='tao_leastsquares_tutorials-chwirut2_1'
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
args='-tao_smonitor -tao_max_it 100 -tao_type pounders -tao_gatol 1.e-5'
diff_args=''
timeoutfactor=1

mpiexec=${PETSCMPIEXEC:-"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/bin/mpiexec"}
diffexec=${PETSCDIFF:-"${petsc_bindir}/petscdiff"}

. "${config_dir}/petsc_harness.sh"

# The diff flags come from script arguments
diff_exe="${diffexec} ${diff_flags} ${diff_args}"
mpiexec="${mpiexec} ${mpiexec_flags}"





nsize_in=${nsize:-3}


for insize in ${nsize_in}; do

   petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args} " chwirut2_1.tmp ${testname}.err "${label}" 
   res=$?

   if test $res = 0; then
      petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/tao/leastsquares/tutorials/output/chwirut2_1.out chwirut2_1.tmp > diff-${testname}-0.out 2> diff-${testname}-0.out || ${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/tao/leastsquares/tutorials/output/chwirut2_1_alt.out chwirut2_1.tmp > diff-${testname}-1.out 2> diff-${testname}-1.out || ${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/tao/leastsquares/tutorials/output/chwirut2_1_alt_2.out chwirut2_1.tmp > diff-${testname}-2.out 2> diff-${testname}-2.out || ${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/tao/leastsquares/tutorials/output/chwirut2_1_alt_3.out chwirut2_1.tmp > diff-${testname}-3.out 2> diff-${testname}-3.out || ${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/tao/leastsquares/tutorials/output/chwirut2_1_alt_4.out chwirut2_1.tmp > diff-${testname}-4.out 2> diff-${testname}-4.out || ${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/tao/leastsquares/tutorials/output/chwirut2_1_alt_5.out chwirut2_1.tmp" diff-${testname}.out diff-${testname}.out diff-${label} ""
   else
      petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
   fi

done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
