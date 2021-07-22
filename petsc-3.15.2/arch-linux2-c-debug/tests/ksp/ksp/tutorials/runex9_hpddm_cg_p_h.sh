#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex9'
testname='runex9_hpddm_cg_p_h'
label='ksp_ksp_tutorials-ex9_hpddm_cg_p_h'
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
args='-unsym 0 -t 2 -pc_type jacobi -ksp_monitor_short -s2_pc_type jacobi -s2_ksp_monitor_short -ksp_rtol 1.e-2 -s2_ksp_rtol 1.e-2 -ksp_type cg -s2_ksp_type hpddm'
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
s2_ksp_hpddm_type_in=${s2_ksp_hpddm_type:-cg bcg bfbcg}


for insize in ${nsize_in}; do
   for is2_ksp_hpddm_type in ${s2_ksp_hpddm_type_in}; do

      petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args}  -s2_ksp_hpddm_type ${is2_ksp_hpddm_type}" ex9_hpddm_cg_p_h.tmp ${testname}.err "${label}+s2_ksp_hpddm_type-${is2_ksp_hpddm_type}" 
      res=$?

      if test $res = 0; then
         petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/ksp/ksp/tutorials/output/ex9_hpddm_cg.out ex9_hpddm_cg_p_h.tmp" diff-${testname}.out diff-${testname}.out diff-${label}+s2_ksp_hpddm_type-${is2_ksp_hpddm_type} ""
      else
         petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
      fi

   done
done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
