#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* sfutils.c */
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

#include "petscsf.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define petscsfsetgraphsection_ PETSCSFSETGRAPHSECTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define petscsfsetgraphsection_ petscsfsetgraphsection
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define petscsfcreatebymatchingindices_ PETSCSFCREATEBYMATCHINGINDICES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define petscsfcreatebymatchingindices_ petscsfcreatebymatchingindices
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
PETSC_EXTERN void  petscsfsetgraphsection_(PetscSF sf,PetscSection localSection,PetscSection globalSection, int *__ierr)
{
*__ierr = PetscSFSetGraphSection(
	(PetscSF)PetscToPointer((sf) ),
	(PetscSection)PetscToPointer((localSection) ),
	(PetscSection)PetscToPointer((globalSection) ));
}
PETSC_EXTERN void  petscsfcreatebymatchingindices_(PetscLayout layout,PetscInt *numRootIndices, PetscInt *rootIndices, PetscInt *rootLocalIndices,PetscInt *rootLocalOffset,PetscInt *numLeafIndices, PetscInt *leafIndices, PetscInt *leafLocalIndices,PetscInt *leafLocalOffset,PetscSF *sfA,PetscSF *sf, int *__ierr)
{
*__ierr = PetscSFCreateByMatchingIndices(
	(PetscLayout)PetscToPointer((layout) ),*numRootIndices,rootIndices,rootLocalIndices,*rootLocalOffset,*numLeafIndices,leafIndices,leafLocalIndices,*leafLocalOffset,sfA,sf);
}
#if defined(__cplusplus)
}
#endif
