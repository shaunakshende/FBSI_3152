#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex6'
testname='runex6_4_ksp_norm_type-natural_pc_type-bjacobi'
label='ksp_ksp_tests-ex6_4_ksp_norm_type-natural_pc_type-bjacobi'
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
args='-ksp_converged_reason -ksp_max_it 20 -ksp_converged_maxits -f ${DATAFILESPATH}/matrices/poisson_2d13p -b_in_f 0 -test_residual -ksp_norm_type natural -pc_type bjacobi'
diff_args=''
timeoutfactor=1

mpiexec=${PETSCMPIEXEC:-"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/bin/mpiexec"}
diffexec=${PETSCDIFF:-"${petsc_bindir}/petscdiff"}

. "${config_dir}/petsc_harness.sh"

# The diff flags come from script arguments
diff_exe="${diffexec} ${diff_flags} ${diff_args}"
mpiexec="${mpiexec} ${mpiexec_flags}"



if ! $force; then
    if test -z "${DATAFILESPATH}"; then
        petsc_report_tapoutput "" "${label}" "SKIP Requires DATAFILESPATH"
        total=1; skip=1
        petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
        exit
    fi
fi


nsize_in=${nsize:-1}
ksp_type_in=${ksp_type:-cg pipecg groppcg}


for insize in ${nsize_in}; do
   for iksp_type in ${ksp_type_in}; do

      petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args}  -ksp_type ${iksp_type}" ex6_4_ksp_norm_type-natural_pc_type-bjacobi.tmp ${testname}.err "${label}+ksp_type-${iksp_type}" 
      res=$?

      if test $res = 0; then
         petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/ksp/ksp/tests/output/ex6_4_ksp_norm_type-natural_pc_type-bjacobi.out ex6_4_ksp_norm_type-natural_pc_type-bjacobi.tmp" diff-${testname}.out diff-${testname}.out diff-${label}+ksp_type-${iksp_type} ""
      else
         petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
      fi

   done
done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
