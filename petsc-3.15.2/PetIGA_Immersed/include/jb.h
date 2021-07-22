
#if !defined(JB_H)
#define JB_H

#include "petiga.h"
#include "../src/petigagrid.h"


//typedef PetscErrorCode (*TSAlphaAdaptFunction)(TS,PetscReal,Vec,Vec,PetscReal*,PetscBool*,void*);


typedef struct {

  PetscReal stage_time;
  PetscReal scale_F;
  PetscReal scale_J;
  PetscReal shift_V;
  PetscReal shift_A;
  Vec       vec_sol_X;
  Vec       vec_sol_V;
  Vec       vec_sol_A;


  Vec       X0,Xa,X1,Xp;
  Vec       V0,Va,V1,Vp,V1_BE, Vnet,Vneta;
  Vec       A0,Aa,A1,Ap,dAs;

  int      (*StageTime)(TS);
  int      (*StageVecs)(TS,Vec);

  TSIFunction2 Function;
  void         *FunCtx;
  TSIJacobian2 Jacobian;
  void         *JacCtx;

  PetscReal Alpha_m;
  PetscReal Alpha_f;
  PetscReal Beta;
  PetscReal Gamma;

  Vec       E;
  TSAlphaAdaptFunction adapt;
  void                 *adaptctx;
  PetscReal            rtol;
  PetscReal            atol;
  PetscReal            rho;
  PetscReal            scale_min;
  PetscReal            scale_max;
  PetscReal            dt_min;
  PetscReal            dt_max;







} TS_Alpha2;

/*
typedef struct {
  PetscReal alpha;        // sufficient decrease parameter
} SNESLineSearch_BT

*/

typedef struct {

//General
  IGA       iga;
  IGA       igaG;
  IGA       igaF;
  PetscInt  Nelt; //Number of solid element
  PetscReal t;      //Width
  PetscReal size;   //Domain size
  PetscReal angle;
  PetscReal g;

  PetscInt  step;
  PetscInt  FreqResults;
  PetscReal TimeRestart;
  PetscInt  StepRestart;
  PetscInt  FreqRestarts;
  PetscReal tolt;
  PetscReal Stmin;

//Fluid
  PetscReal Tcte,acte,bcte,Rcte;
  PetscReal lambda;
  PetscReal Cv,visc,Kcond;
  PetscReal T_wall,T_ini,rho_ini;

//Mesh
  PetscReal lambda_mesh,mu_mesh;

//Solid
  PetscReal kappa_str,mu_str,rho_str;
  PetscReal pri;

//  void (*model) (IGAPoint pnt, const PetscScalar *U, PetscScalar (*P)[1], PetscScalar (*S)[1], PetscScalar (*D)[1][1][1], void *ctx);
//  void (*model) (IGAPoint pnt, const PetscScalar *U, PetscScalar (*P)[2], PetscScalar (*S)[2], PetscScalar (*D)[2][2][2], void *ctx);
  void (*model) (IGAPoint pnt, const PetscScalar *U, PetscScalar (*P)[2], PetscScalar (*S)[2], PetscScalar (*D)[2][2][2], void *ctx, PetscScalar (*Finv)[2], PetscReal *J);
//  void (*model) (IGAPoint pnt, const PetscScalar *U, PetscScalar (*P)[3], PetscScalar (*S)[3], PetscScalar (*D)[3][3][3], void *ctx);


  // problem parameters
  PetscReal L0,h;
  PetscReal Ca,alpha,theta,Re;
  // bubble centers
  PetscReal C1[2],R1;
  PetscReal C2[2],R2;
  PetscReal C3[2],R3;
  PetscReal C4[2],R4;
  PetscReal C5[2],R5;
  PetscReal C6[2],R6;
  PetscReal C7[2],R7;



} AppCtx;


/*
#undef  __FUNCT__
#define __FUNCT__ "IGA_NewGridIO"
static PetscErrorCode IGA_NewGridIO(IGA iga,PetscInt bs,IGA_Grid *grid)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidPointer(grid,3);
  {
    MPI_Comm comm    = ((PetscObject)iga)->comm;
    PetscInt dim     = iga->dim;
    PetscInt *sizes  = iga->geom_sizes;
    PetscInt *lstart = iga->geom_lstart;
    PetscInt *lwidth = iga->geom_lwidth;
    PetscInt *gstart = iga->geom_gstart;
    PetscInt *gwidth = iga->geom_gwidth;
    ierr = IGA_Grid_Create(comm,grid);CHKERRQ(ierr);
    ierr = IGA_Grid_Init(*grid,dim,bs,sizes,lstart,lwidth,gstart,gwidth);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

*/







#endif/*JB_H*/
