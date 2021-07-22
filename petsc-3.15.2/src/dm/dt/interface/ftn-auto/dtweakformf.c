#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* dtweakform.c */
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

#include "petscds.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define petscweakformgetnumfields_ PETSCWEAKFORMGETNUMFIELDS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define petscweakformgetnumfields_ petscweakformgetnumfields
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define petscweakformsetnumfields_ PETSCWEAKFORMSETNUMFIELDS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define petscweakformsetnumfields_ petscweakformsetnumfields
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define petscweakformdestroy_ PETSCWEAKFORMDESTROY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define petscweakformdestroy_ petscweakformdestroy
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define petscweakformcreate_ PETSCWEAKFORMCREATE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define petscweakformcreate_ petscweakformcreate
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
PETSC_EXTERN void  petscweakformgetnumfields_(PetscWeakForm *wf,PetscInt *Nf, int *__ierr)
{
*__ierr = PetscWeakFormGetNumFields(*wf,Nf);
}
PETSC_EXTERN void  petscweakformsetnumfields_(PetscWeakForm *wf,PetscInt *Nf, int *__ierr)
{
*__ierr = PetscWeakFormSetNumFields(*wf,*Nf);
}
PETSC_EXTERN void  petscweakformdestroy_(PetscWeakForm *wf, int *__ierr)
{
*__ierr = PetscWeakFormDestroy(
	(PetscWeakForm* )PetscToPointer((wf) ));
}
PETSC_EXTERN void  petscweakformcreate_(MPI_Fint * comm,PetscWeakForm *wf, int *__ierr)
{
*__ierr = PetscWeakFormCreate(
	MPI_Comm_f2c(*(comm)),
	(PetscWeakForm* )PetscToPointer((wf) ));
}
#if defined(__cplusplus)
}
#endif
