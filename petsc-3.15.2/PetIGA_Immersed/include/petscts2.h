#if !defined(__PETSCTS2_H)
#define __PETSCTS2_H
#include <petscts.h>


typedef PetscErrorCode (*TSAlphaAdaptFunction)(TS,PetscReal,Vec,Vec,PetscReal*,PetscBool*,void*);

typedef PetscErrorCode (*TSIFunction2)(TS,PetscReal,Vec,Vec,Vec,Vec,void*);
typedef PetscErrorCode (*TSIJacobian2)(TS,PetscReal,Vec,Vec,Vec,PetscReal,PetscReal,Mat*,Mat*,MatStructure*,void*);
PETSC_EXTERN PetscErrorCode TSSetIFunction2(TS,Vec,TSIFunction2,void*);
PETSC_EXTERN PetscErrorCode TSSetIJacobian2(TS,Mat,Mat,TSIJacobian2,void*);
PETSC_EXTERN PetscErrorCode TSComputeIFunction2(TS,PetscReal,Vec,Vec,Vec,Vec,PetscBool);
PETSC_EXTERN PetscErrorCode TSComputeIJacobian2(TS,PetscReal,Vec,Vec,Vec,PetscReal,PetscReal,Mat*,Mat*,MatStructure*,PetscBool);
PETSC_EXTERN PetscErrorCode TSSetSolution2(TS,Vec,Vec);
PETSC_EXTERN PetscErrorCode TSGetSolution2(TS,Vec*,Vec*);
PETSC_EXTERN PetscErrorCode TSSolve2(TS,Vec,Vec);

#define TSALPHA2 "alpha2"
PETSC_EXTERN PetscErrorCode TSAlpha2SetRadius(TS,PetscReal);
PETSC_EXTERN PetscErrorCode TSAlpha2SetParams(TS,PetscReal,PetscReal,PetscReal,PetscReal);
PETSC_EXTERN PetscErrorCode TSAlpha2GetParams(TS,PetscReal*,PetscReal*,PetscReal*,PetscReal*);

PETSC_EXTERN PetscErrorCode TSAlpha2AdaptDefault_JB(TS,PetscReal,Vec,Vec, PetscReal*,PetscBool *,void*);
PETSC_EXTERN PetscErrorCode TSAlpha2SetAdapt_JB(TS,TSAlphaAdaptFunction,void*);


#endif/*__PETSCTS2_H*/
