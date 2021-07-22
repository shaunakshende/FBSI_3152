#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex6'
testname='runex6_3_skip_hpddm'
label='ksp_ksp_tests-ex6_3_skip_hpddm'
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
args='-f ${wPETSC_DIR}/share/petsc/datafiles/matrices/spd-real-int32-float64 -pc_type none -ksp_max_it 8 -ksp_error_if_not_converged -ksp_convergence_test skip -ksp_converged_reason -test_residual -ksp_type hpddm'
diff_args=''
timeoutfactor=1

mpiexec=${PETSCMPIEXEC:-"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/bin/mpiexec"}
diffexec=${PETSCDIFF:-"${petsc_bindir}/petscdiff"}

. "${config_dir}/petsc_harness.sh"

# The diff flags come from script arguments
diff_exe="${diffexec} ${diff_flags} ${diff_args}"
mpiexec="${mpiexec} ${mpiexec_flags}"



if ! $force; then
    petsc_report_tapoutput "" "${label}" "SKIP PETSC_HAVE_HPDDM requirement not met"
    total=1; skip=1
    petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
    exit
fi


nsize_in=${nsize:-1}
ksp_hpddm_type_in=${ksp_hpddm_type:-cg gmres bgmres bcg bfbcg gcrodr bgcrodr}


for insize in ${nsize_in}; do
   for iksp_hpddm_type in ${ksp_hpddm_type_in}; do

      petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args}  -ksp_hpddm_type ${iksp_hpddm_type}" ex6_3_skip_hpddm.tmp ${testname}.err "${label}+ksp_hpddm_type-${iksp_hpddm_type}" 
      res=$?

      if test $res = 0; then
         petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/ksp/ksp/tests/output/ex6_skip.out ex6_3_skip_hpddm.tmp" diff-${testname}.out diff-${testname}.out diff-${label}+ksp_hpddm_type-${iksp_hpddm_type} ""
      else
         petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
      fi

   done
done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
