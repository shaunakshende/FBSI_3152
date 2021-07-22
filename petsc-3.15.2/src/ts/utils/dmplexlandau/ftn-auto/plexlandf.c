#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* plexland.c */
/* Fortran interface file */

/*
* This file was generated automatically by bfort from the C source
* file.  
 */

#ifdef PETSC_USE_POINTER_CONVERSION
#if defined(__cplusplus)
extern "C" { 
#endif 
extern void *PetscToPointer(void*);
extern int PetscFromPointer(void *);
extern void PetscRmPointer(void*);
#if defined(__cplusplus)
} 
#endif 

#else

#define PetscToPointer(a) (*(PetscFortranAddr *)(a))
#define PetscFromPointer(a) (PetscFortranAddr)(a)
#define PetscRmPointer(a)
#endif

#include "petscdmplex.h"
#include "petsclandau.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define landauaddmaxwellians_ LANDAUADDMAXWELLIANS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define landauaddmaxwellians_ landauaddmaxwellians
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define landaudestroyvelocityspace_ LANDAUDESTROYVELOCITYSPACE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define landaudestroyvelocityspace_ landaudestroyvelocityspace
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define landauprintnorms_ LANDAUPRINTNORMS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define landauprintnorms_ landauprintnorms
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define landaucreatecoloring_ LANDAUCREATECOLORING
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define landaucreatecoloring_ landaucreatecoloring
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define landaucreatemassmatrix_ LANDAUCREATEMASSMATRIX
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define landaucreatemassmatrix_ landaucreatemassmatrix
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define landauifunction_ LANDAUIFUNCTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define landauifunction_ landauifunction
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define landauijacobian_ LANDAUIJACOBIAN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define landauijacobian_ landauijacobian
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
PETSC_EXTERN void  landauaddmaxwellians_(DM dm,Vec X,PetscReal *time,PetscReal temps[],PetscReal ns[],void*actx, int *__ierr)
{
*__ierr = LandauAddMaxwellians(
	(DM)PetscToPointer((dm) ),
	(Vec)PetscToPointer((X) ),*time,temps,ns,actx);
}
PETSC_EXTERN void  landaudestroyvelocityspace_(DM *dm, int *__ierr)
{
*__ierr = LandauDestroyVelocitySpace(dm);
}
PETSC_EXTERN void  landauprintnorms_(Vec X,PetscInt *stepi, int *__ierr)
{
*__ierr = LandauPrintNorms(
	(Vec)PetscToPointer((X) ),*stepi);
}
PETSC_EXTERN void  landaucreatecoloring_(Mat JacP,DM plex,PetscContainer *container, int *__ierr)
{
*__ierr = LandauCreateColoring(
	(Mat)PetscToPointer((JacP) ),
	(DM)PetscToPointer((plex) ),container);
}
PETSC_EXTERN void  landaucreatemassmatrix_(DM dm,Mat *Amat, int *__ierr)
{
*__ierr = LandauCreateMassMatrix(
	(DM)PetscToPointer((dm) ),Amat);
}
PETSC_EXTERN void  landauifunction_(TS ts,PetscReal *time_dummy,Vec X,Vec X_t,Vec F,void*actx, int *__ierr)
{
*__ierr = LandauIFunction(
	(TS)PetscToPointer((ts) ),*time_dummy,
	(Vec)PetscToPointer((X) ),
	(Vec)PetscToPointer((X_t) ),
	(Vec)PetscToPointer((F) ),actx);
}
PETSC_EXTERN void  landauijacobian_(TS ts,PetscReal *time_dummy,Vec X,Vec U_tdummy,PetscReal *shift,Mat Amat,Mat Pmat,void*actx, int *__ierr)
{
*__ierr = LandauIJacobian(
	(TS)PetscToPointer((ts) ),*time_dummy,
	(Vec)PetscToPointer((X) ),
	(Vec)PetscToPointer((U_tdummy) ),*shift,
	(Mat)PetscToPointer((Amat) ),
	(Mat)PetscToPointer((Pmat) ),actx);
}
#if defined(__cplusplus)
}
#endif
