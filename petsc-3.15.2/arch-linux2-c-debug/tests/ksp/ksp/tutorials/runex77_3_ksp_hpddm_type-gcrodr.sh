#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex77'
testname='runex77_3_ksp_hpddm_type-gcrodr'
label='ksp_ksp_tutorials-ex77_3_ksp_hpddm_type-gcrodr'
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
args='-ksp_converged_reason -ksp_max_it 500 -f ${DATAFILESPATH}/matrices/hpddm/GCRODR/A_400.dat -ksp_type hpddm -ksp_hpddm_recycle 5 -ksp_hpddm_type gcrodr'
diff_args=''
timeoutfactor=1

mpiexec=${PETSCMPIEXEC:-"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/bin/mpiexec"}
diffexec=${PETSCDIFF:-"${petsc_bindir}/petscdiff"}

. "${config_dir}/petsc_harness.sh"

# The diff flags come from script arguments
diff_exe="${diffexec} ${diff_flags} ${diff_args}"
mpiexec="${mpiexec} ${mpiexec_flags}"



if ! $force; then
    petsc_report_tapoutput "" "${label}" "SKIP Requires DATAFILESPATH, PETSC_HAVE_HPDDM requirement not met"
    total=1; skip=1
    petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
    exit
fi


nsize_in=${nsize:-2}
mat_type_in=${mat_type:-aij sbaij}


for insize in ${nsize_in}; do
   for imat_type in ${mat_type_in}; do

      petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args}  -mat_type ${imat_type}" ex77_3_ksp_hpddm_type-gcrodr.tmp ${testname}.err "${label}+mat_type-${imat_type}" 
      res=$?

      if test $res = 0; then
         petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/ksp/ksp/tutorials/output/ex77_3_ksp_hpddm_type-gcrodr.out ex77_3_ksp_hpddm_type-gcrodr.tmp" diff-${testname}.out diff-${testname}.out diff-${label}+mat_type-${imat_type} ""
      else
         petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
      fi

   done
done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
