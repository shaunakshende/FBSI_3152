#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex33'
testname='runex33_parmetis_nsize-1_nparts-2'
label='dm_partitioner_tests-ex33_parmetis_nsize-1_nparts-2'
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
args='-petscpartitioner_type parmetis -petscpartitioner_view -petscpartitioner_view_graph -nparts 2'
diff_args=''
timeoutfactor=1

mpiexec=${PETSCMPIEXEC:-"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/bin/mpiexec"}
diffexec=${PETSCDIFF:-"${petsc_bindir}/petscdiff"}

. "${config_dir}/petsc_harness.sh"

# The diff flags come from script arguments
diff_exe="${diffexec} ${diff_flags} ${diff_args}"
mpiexec="${mpiexec} ${mpiexec_flags}"



if ! $force; then
    petsc_report_tapoutput "" "${label}" "SKIP PETSC_HAVE_PARMETIS requirement not met"
    total=1; skip=1
    petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
    exit
fi


nsize_in=${nsize:-1}
pwgts_in=${pwgts:-false true}
vwgts_in=${vwgts:-false true}


for insize in ${nsize_in}; do
   for ipwgts in ${pwgts_in}; do
      for ivwgts in ${vwgts_in}; do

         petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args}  -pwgts ${ipwgts} -vwgts ${ivwgts}" ex33_parmetis_nsize-1_nparts-2.tmp ${testname}.err "${label}+pwgts-${ipwgts}_vwgts-${ivwgts}" 
         res=$?

         if test $res = 0; then
            petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/dm/partitioner/tests/output/ex33_parmetis_nsize-1_nparts-2.out ex33_parmetis_nsize-1_nparts-2.tmp" diff-${testname}.out diff-${testname}.out diff-${label}+pwgts-${ipwgts}_vwgts-${ivwgts} ""
         else
            petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
         fi

      done
   done
done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 