
#requiresdefine   'PETSC_HAVE_FORTRAN'
ALL: lib

CFLAGS   =
FFLAGS   =
SOURCEC  = sortsof.c mpitsf.c sortif.c sortipf.c sortdf.c psplitf.c mpiuf.c memcf.c psleepf.c
SOURCEF  =
SOURCEH  = 
DIRS     =
LIBBASE	  = libpetscsys
LOCDIR	  = src/sys/utils/ftn-auto/
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test
