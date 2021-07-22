#if !defined(GLOBAL_H)
#define GLOBAL_H

/*
#include <petscconf.h>
#undef  PETSC_STATIC_INLINE
#define PETSC_STATIC_INLINE static __inline
*/
#include <petiga.h>
#include <petsc.h>
#include <petsc-private/petscimpl.h>
#include <petscts2.h>



/* ---------------------------------------------------------------- */


//typedef struct _n_IGAAxis     *IGAAxis;
typedef struct _n_RKNode     *RKNode;

/* ---------------------------------------------------------------- */

//struct _n_IGAAxis {
//  PetscInt refct;
//  /**/
//  PetscInt   p; /* polynomial order    */
//  PetscInt   m; /* last knot index     */
//  PetscReal *U; /* knot vector         */
//  /**/
//  PetscBool  periodic; /* periodicity  */
//  PetscInt   nnp,nel;  /* bases, spans */
//  PetscInt   *span;    /* span indices */
//};
//PETSC_EXTERN PetscErrorCode IGAAxisCreate(IGAAxis *axis);
//PETSC_EXTERN PetscErrorCode IGAAxisDestroy(IGAAxis *axis);
//PETSC_EXTERN PetscErrorCode IGAAxisReset(IGAAxis axis);
//PETSC_EXTERN PetscErrorCode IGAAxisReference(IGAAxis axis);
//PETSC_EXTERN PetscErrorCode IGAAxisCopy(IGAAxis base,IGAAxis axis);
//PETSC_EXTERN PetscErrorCode IGAAxisDuplicate(IGAAxis base,IGAAxis *axis);
//PETSC_EXTERN PetscErrorCode IGAAxisSetPeriodic(IGAAxis axis,PetscBool periodic);
//PETSC_EXTERN PetscErrorCode IGAAxisGetPeriodic(IGAAxis axis,PetscBool *periodic);
//PETSC_EXTERN PetscErrorCode IGAAxisSetDegree(IGAAxis axis,PetscInt p);
//PETSC_EXTERN PetscErrorCode IGAAxisGetDegree(IGAAxis axis,PetscInt *p);
//PETSC_EXTERN PetscErrorCode IGAAxisSetKnots(IGAAxis axis,PetscInt m,const PetscReal U[]);
//PETSC_EXTERN PetscErrorCode IGAAxisGetKnots(IGAAxis axis,PetscInt *m,PetscReal *U[]);
//PETSC_EXTERN PetscErrorCode IGAAxisGetLimits(IGAAxis axis,PetscReal *Ui,PetscReal *Uf);
//PETSC_EXTERN PetscErrorCode IGAAxisGetSizes(IGAAxis axis,PetscInt *nel,PetscInt *nnp);
//PETSC_EXTERN PetscErrorCode IGAAxisGetSpans(IGAAxis axis,PetscInt *nel,PetscInt *spans[]);
//PETSC_EXTERN PetscErrorCode IGAAxisInit(IGAAxis axis,PetscInt p,PetscInt m,const PetscReal U[]);
//PETSC_EXTERN PetscErrorCode IGAAxisInitBreaks(IGAAxis axis,PetscInt nu,const PetscReal u[],PetscInt C);
//PETSC_EXTERN PetscErrorCode IGAAxisInitUniform(IGAAxis axis,PetscInt N,PetscReal Ui,PetscReal Uf,PetscInt C);
//PETSC_EXTERN PetscErrorCode IGAAxisSetUp(IGAAxis axis);



struct _n_RKNode {

  unsigned long int  nodeID;
  double         support[3];
  double         supportInverse[3];
  double         nodalVolume;
  double         integrationNodalVolume;
  PetscInt         nodalMass;

  double         nodalForce[3];
  double         acc0[3];
  long int       maskingArray[3];
//  vector<FacePoint>         gaussIntegrationPoints;
//  vector <vector <double> > smoothedGradient;

  double         effectiveModulus;
  double         currentDeformationGradient[9];
  double         incrementalDeformationGradient[9];
//  vector < vector<double> > smoothedSpatialGradient;

  double         determinantDeformationGradient;
  double         totalGeneralizedDisplacement[3];
  double         incrementalGeneralizedDisplacement[3];

  double         incrementalGeneralizedVelocity[3];
  double         totalGeneralizedVelocity[3];

  PetscScalar  referenceCoord[3];

  int            isDirichletBoundary;
  double         prescribedIncrementalDisplacements[3];
  double         prescribedTotalDisplacements[3];


};
PETSC_EXTERN PetscErrorCode preProcess();


//  vector <double>    computeShapeFunc(const Parameters &parameters) ;
//PETSC_EXTERN PetscErrorCode  computeSmoothedGradient(const Parameters &parameters);
//  void           computeSmoothedSpatialGradient(const Parameters &parameters);
//  void           computeInternalForces(Parameters & __restrict__ parameters,vector <RKNode*> & __restrict__ tableOfNodes);
//  void           getStresses(const Parameters &parameters);
//  void           computeDeformationGradients( const double *__restrict__ Finit, bool skipDefGrad);
//  void           createSNNIIntegrationPoints();











/* ---------------------------------------------------------------- */

#if defined(PETSC_USE_DEBUG)
#define IGACheckSetUp(iga,arg) do {                                      \
    if (PetscUnlikely(!(iga)->setup))                                    \
      SETERRQ3(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,                 \
               "Must call IGASetUp() on argument %D \"%s\" before %s()", \
               (arg),#iga,PETSC_FUNCTION_NAME);                          \
    } while (0)
#define IGACheckSetUpStage(iga,arg,stg) do {                             \
    if (PetscUnlikely((iga)->setupstage<(stg))) IGACheckSetUp(iga,arg);  \
    } while (0)
#else
#define IGACheckSetUp(iga,arg)          do {} while (0)
#define IGACheckSetUpStage(iga,arg,stg) do {} while (0)
#endif
#define IGACheckSetUpStage1(iga,arg) IGACheckSetUpStage(iga,arg,1)
#define IGACheckSetUpStage2(iga,arg) IGACheckSetUpStage(iga,arg,2)

#if defined(PETSC_USE_DEBUG)
#define IGACheckFormOp(iga,arg,FormOp) do {             \
    if (!iga->form->ops->FormOp)                        \
      SETERRQ4(((PetscObject)iga)->comm,PETSC_ERR_USER, \
               "Must call IGASetForm%s() "              \
               "on argument %D \"%s\" before %s()",     \
               #FormOp,(arg),#iga,PETSC_FUNCTION_NAME); \
    } while (0)
#else
#define IGACheckFormOp(iga,arg,FormOp) do {} while (0)
#endif

/* ---------------------------------------------------------------- */

#ifndef PetscValidRealPointer
#define PetscValidRealPointer PetscValidDoublePointer
#endif

#ifndef PetscMalloc1
#define PetscMalloc1(m1,r1) \
  PetscMalloc((m1)*sizeof(**(r1)),r1)
#endif

#ifndef PetscCalloc1
#define PetscCalloc1(m1,r1) \
  (PetscMalloc1((m1),PetscInt,r1) || PetscMemzero(*(r1),(m1)*sizeof(**(r1))))
#endif

/* ---------------------------------------------------------------- */

#if PETSC_VERSION_(3,3,0)

#undef  PETSC_VERSION_LT
#define PETSC_VERSION_LT(MAJOR,MINOR,SUBMINOR)    \
  (PETSC_VERSION_RELEASE == 1 &&                  \
   (PETSC_VERSION_MAJOR < (MAJOR) ||              \
    (PETSC_VERSION_MAJOR == (MAJOR) &&            \
     (PETSC_VERSION_MINOR < (MINOR) ||            \
      (PETSC_VERSION_MINOR == (MINOR) &&          \
       (PETSC_VERSION_SUBMINOR < (SUBMINOR)))))))

#undef  PETSC_VERSION_LE
#define PETSC_VERSION_LE(MAJOR,MINOR,SUBMINOR) \
  (PETSC_VERSION_LT(MAJOR,MINOR,SUBMINOR) ||   \
   PETSC_VERSION_(MAJOR,MINOR,SUBMINOR))

#undef  PETSC_VERSION_GE
#define PETSC_VERSION_GE(MAJOR,MINOR,SUBMINOR) \
  (0 == PETSC_VERSION_LT(MAJOR,MINOR,SUBMINOR))

#undef  PETSC_VERSION_GT
#define PETSC_VERSION_GT(MAJOR,MINOR,SUBMINOR) \
  (0 == PETSC_VERSION_LE(MAJOR,MINOR,SUBMINOR))

#endif

/* ---------------------------------------------------------------- */

#endif/*PETIGA_H*/
