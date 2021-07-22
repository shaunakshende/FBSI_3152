#!/usr/bin/env bash
# This script was created by gmakegentest.py



# PATH for DLLs on windows
PATH="$PATH":"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/lib"

exec='../ex99'
testname='runex99_msh-vtx'
label='dm_impls_plex_tests-ex99_msh-vtx'
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
args='-dir ${wPETSC_DIR}/share/petsc/datafiles/meshes -order 1 -dm_view ::ascii_info_detail -dm_plex_check_all -dm_plex_gmsh_highorder false -dm_plex_gmsh_use_marker true -dm_plex_gmsh_spacedim 3 -msh vtx'
diff_args=''
timeoutfactor=1

mpiexec=${PETSCMPIEXEC:-"/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/bin/mpiexec"}
diffexec=${PETSCDIFF:-"${petsc_bindir}/petscdiff"}

. "${config_dir}/petsc_harness.sh"

# The diff flags come from script arguments
diff_exe="${diffexec} ${diff_flags} ${diff_args}"
mpiexec="${mpiexec} ${mpiexec_flags}"



if ! $force; then
    petsc_report_tapoutput "" "${label}" "SKIP Required: define(PETSC_GMSH_EXE)"
    total=1; skip=1
    petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
    exit
fi


nsize_in=${nsize:-1}
fmt_in=${fmt:-msh22 msh40 msh41}
bin_in=${bin:-0 1}


for insize in ${nsize_in}; do
   for ifmt in ${fmt_in}; do
      for ibin in ${bin_in}; do

         petsc_testrun "${mpiexec} -n ${insize} ${exec} ${args}  -fmt ${ifmt} -bin ${ibin}" ex99_msh-vtx.tmp ${testname}.err "${label}+fmt-${ifmt}_bin-${ibin}" 
         res=$?

         if test $res = 0; then
            petsc_testrun "${diff_exe} /home/shaunak/Desktop/petsc-3.15.2/src/dm/impls/plex/tests/output/ex99_msh-vtx.out ex99_msh-vtx.tmp" diff-${testname}.out diff-${testname}.out diff-${label}+fmt-${ifmt}_bin-${ibin} ""
         else
            petsc_report_tapoutput "" ${label} "SKIP Command failed so no diff"
         fi

      done
   done
done

petsc_testend "/home/shaunak/Desktop/petsc-3.15.2/arch-linux2-c-debug/tests" 
