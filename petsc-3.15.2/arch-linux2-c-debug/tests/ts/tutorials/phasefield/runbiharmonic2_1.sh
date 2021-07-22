#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../biharmonic2'
testname='runbiharmonic2_1'
label='ts_tutorials_phasefield-biharmonic2_1'
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
args='-ts_monitor -snes_monitor -pc_type lu -snes_converged_reason -ts_type beuler -da_refine 5 -ts_dt 9.53674e-9 -ts_max_steps 50'
diff_args=''
timeoutfactor=1

mpiexec=${PETSCMPIEXEC:-"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/bin/mpiexec"}
diffexec=${PETSCDIFF:-"${petsc_bindir}/petscdiff"}

. "${config_dir}/petsc_harness.sh"

# The diff flags come from script arguments
diff_exe="${diffexec} ${diff_flags} ${diff_args}"
mpiexec="${mpiexec} ${mpiexec_flags}"



if ! $force; then
    petsc_report_tapoutput "" "${label}" "SKIP PETSC_HAVE_X requirement not met"
    total=1; skip=1
    petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
    exit
fi


nsize_in=${nsize:-1}


for insize in ${nsize_in}; do

   petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args} " biharmonic2_1.tmp ${testname}.err "${label}" 
   res=$?

   if test $res = 0; then
      petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/ts/tutorials/phasefield/output/biharmonic2_1.out biharmonic2_1.tmp > diff-${testname}-0.out 2> diff-${testname}-0.out || ${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/ts/tutorials/phasefield/output/biharmonic2_1_alt.out biharmonic2_1.tmp > diff-${testname}-1.out 2> diff-${testname}-1.out || ${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/ts/tutorials/phasefield/output/biharmonic2_1_alt_2.out biharmonic2_1.tmp > diff-${testname}-2.out 2> diff-${testname}-2.out || ${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/ts/tutorials/phasefield/output/biharmonic2_1_alt_3.out biharmonic2_1.tmp > diff-${testname}-3.out 2> diff-${testname}-3.out || ${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/ts/tutorials/phasefield/output/biharmonic2_1_alt_4.out biharmonic2_1.tmp" diff-${testname}.out diff-${testname}.out diff-${label} ""
   else
      petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
   fi

done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
