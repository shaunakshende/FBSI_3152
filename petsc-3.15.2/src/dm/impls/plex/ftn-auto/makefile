
#requiresdefine   'PETSC_HAVE_FORTRAN'
ALL: lib
CPPFLAGS = ${NETCFD_INCLUDE} ${EXODUSII_INCLUDE}
CFLAGS   =
FFLAGS   =
SOURCEC  = plexreorderf.c plexpartitionf.c plexegadsf.c plexprojectf.c plexcheckinterfacef.c plexpointf.c plexexodusiif.c plexnaturalf.c plexdistributef.c plexorientf.c plexrefinef.c plexcgnsf.c plexcreatef.c plexfemf.c plexf.c plexgeometryf.c plextreef.c plexfvmf.c plexgmshf.c plexsubmeshf.c plexindicesf.c plexinterpolatef.c
SOURCEF  =
SOURCEH  = 
DIRS     =
LIBBASE  = libpetscdm
LOCDIR   = src/dm/impls/plex/ftn-auto/
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test
