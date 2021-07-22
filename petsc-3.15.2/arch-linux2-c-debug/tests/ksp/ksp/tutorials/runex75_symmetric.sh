#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex75'
testname='runex75_symmetric'
label='ksp_ksp_tutorials-ex75_symmetric'
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
args='-nmat 3 -pc_type jacobi -ksp_converged_reason -ksp_type hpddm -ksp_max_it 1000 -ksp_gmres_restart 40 -ksp_atol 1e-11 -ksp_hpddm_type bgcrodr -ksp_hpddm_recycle 20 -load_dir ${DATAFILESPATH}/matrices/hpddm/GCRODR -ksp_hpddm_recycle_symmetric true'
diff_args=''
timeoutfactor=1

mpiexec=${PETSCMPIEXEC:-"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/bin/mpiexec"}
diffexec=${PETSCDIFF:-"${petsc_bindir}/petscdiff"}

. "${config_dir}/petsc_harness.sh"

# The diff flags come from script arguments
diff_exe="${diffexec} ${diff_flags} ${diff_args}"
mpiexec="${mpiexec} ${mpiexec_flags}"



if ! $force; then
    petsc_report_tapoutput "" "${label}" "SKIP PETSC_HAVE_HPDDM requirement not met, Requires DATAFILESPATH, PETSC_HAVE_SLEPC requirement not met"
    total=1; skip=1
    petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
    exit
fi


nsize_in=${nsize:-2}
reset_in=${reset:-false true}


for insize in ${nsize_in}; do
   for ireset in ${reset_in}; do

      petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args}  -reset ${ireset}" ex75_symmetric.tmp ${testname}.err "${label}+reset-${ireset}" 
      res=$?

      if test $res = 0; then
         petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/ksp/ksp/tutorials/output/ex75_symmetric.out ex75_symmetric.tmp" diff-${testname}.out diff-${testname}.out diff-${label}+reset-${ireset} ""
      else
         petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
      fi

   done
done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 