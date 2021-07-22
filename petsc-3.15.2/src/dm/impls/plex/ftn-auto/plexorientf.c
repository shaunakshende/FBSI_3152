#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* plexorient.c */
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
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmplexcompareorientations_ DMPLEXCOMPAREORIENTATIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmplexcompareorientations_ dmplexcompareorientations
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmplexorientcell_ DMPLEXORIENTCELL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmplexorientcell_ dmplexorientcell
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmplexreversecell_ DMPLEXREVERSECELL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmplexreversecell_ dmplexreversecell
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmplexorient_ DMPLEXORIENT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmplexorient_ dmplexorient
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
PETSC_EXTERN void  dmplexcompareorientations_(DM dm,PetscInt *p,PetscInt *mainConeSize, PetscInt mainCone[],PetscInt *start,PetscBool *reverse, int *__ierr)
{
*__ierr = DMPlexCompareOrientations(
	(DM)PetscToPointer((dm) ),*p,*mainConeSize,mainCone,start,reverse);
}
PETSC_EXTERN void  dmplexorientcell_(DM dm,PetscInt *p,PetscInt *mainConeSize, PetscInt mainCone[], int *__ierr)
{
*__ierr = DMPlexOrientCell(
	(DM)PetscToPointer((dm) ),*p,*mainConeSize,mainCone);
}
PETSC_EXTERN void  dmplexreversecell_(DM dm,PetscInt *cell, int *__ierr)
{
*__ierr = DMPlexReverseCell(
	(DM)PetscToPointer((dm) ),*cell);
}
PETSC_EXTERN void  dmplexorient_(DM dm, int *__ierr)
{
*__ierr = DMPlexOrient(
	(DM)PetscToPointer((dm) ));
}
#if defined(__cplusplus)
}
#endif
