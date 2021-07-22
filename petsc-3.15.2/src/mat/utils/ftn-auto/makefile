
#requiresdefine   'PETSC_HAVE_FORTRAN'
ALL: lib

CFLAGS   =
FFLAGS   =
SOURCEC  = zerodiagf.c gcreatef.c multequalf.c getcolvf.c axpyf.c bandwidthf.c
SOURCEF  =
SOURCEH  = 
DIRS     =
LIBBASE  = libpetscmat
LOCDIR   = src/mat/utils/ftn-auto/
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test
