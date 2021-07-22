#include "petiga.h"
#include "petscts2.h"
#include <petsc-private/tsimpl.h>
#include "jb.h"


PETSC_STATIC_INLINE
PetscBool IGAElementNextUserIFunction(IGAElement element,IGAUserIFunction2 *fun,void **ctx)
{
  IGAUserOps ops;
  while (IGAElementNextUserOps(element,&ops) && !ops->IFunction2);
  if (!ops) return PETSC_FALSE;
  *fun = ops->IFunction2;
  *ctx = ops->IFunCtx;
  return PETSC_TRUE;
}

PETSC_STATIC_INLINE
PetscBool IGAElementNextUserIJacobian(IGAElement element,IGAUserIJacobian2 *jac,void **ctx)
{
  IGAUserOps ops;
  while (IGAElementNextUserOps(element,&ops) && !ops->IJacobian2);
  if (!ops) return PETSC_FALSE;
  *jac = ops->IJacobian2;
  *ctx = ops->IJacCtx;
  return PETSC_TRUE;
}

#undef  __FUNCT__
#define __FUNCT__ "IGAComputeIFunction2"
PetscErrorCode IGAComputeIFunction2(IGA iga,PetscReal dt,
                                    PetscReal a,Vec vecA,
                                    PetscReal v,Vec vecV,
                                    PetscReal t,Vec vecU,
                                    Vec vecF, TS ts)     //TS ts
{
  Vec               localA;
  Vec               localV;
  Vec               localU;
  Vec				localVnet;
  const PetscScalar *arrayA;
  const PetscScalar *arrayV;
  const PetscScalar *arrayU;
  const PetscScalar *arrayVnet;
  IGAElement        element;
  IGAPoint          point;
  IGAUserIFunction2    IFunction;
  void              *ctx;
  PetscScalar       *A,*V,*U,*F,*R, *Vn;
  PetscErrorCode    ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidHeaderSpecific(vecA,VEC_CLASSID,4);
  PetscValidHeaderSpecific(vecV,VEC_CLASSID,6);
  PetscValidHeaderSpecific(vecU,VEC_CLASSID,8);
  PetscValidHeaderSpecific(vecF,VEC_CLASSID,9);
  IGACheckSetUp(iga,1);
  IGACheckUserOp(iga,1,IFunction2);
  TS_Alpha2     *th    = (TS_Alpha2*)ts->data;


  if (iga->cruz == 1){           //(Xag - X0g) Solve for a DeltaX

	  PetscInt i,j,Nloc;
	  PetscInt dof = iga->dof;
	  PetscInt dim = iga->dim;
	  PetscScalar *arrayUa, *arrayX0g;
	  ierr  = VecGetArray(vecU,&arrayUa);CHKERRQ(ierr);
	  ierr  = VecGetArray(th->X0,&arrayX0g);CHKERRQ(ierr);
	  ierr  = VecGetLocalSize(th->X0,&Nloc);CHKERRQ(ierr);

//	  	 PetscInt k
//	     for(i=0;i<15;i++) {
	  //	   	  for(k=0;k<2;k++) {
//	   		    PetscPrintf(PETSC_COMM_WORLD,"i: %d, arrayX1: %e, arrayV1: %e, arrayA1: %e \n",i, arrayX1[i*2+k],arrayV1[i*2+k],arrayA1[i*2+k]);
//	   	  }
//	     }

	     for(i=0;i<Nloc/dof;i++) {
	      for(j=0;j<dim;j++)  arrayUa[dof*i+j] -= arrayX0g[dof*i+j];
	     }

	   VecRestoreArray(vecU,&arrayUa);CHKERRQ(ierr);
	   VecRestoreArray(th->X0,&arrayX0g);CHKERRQ(ierr);
  }






  /* Clear global vector F */
  ierr = VecZeroEntries(vecF);CHKERRQ(ierr);

  /* Get local vectors A,V,U and arrays */
  ierr = IGAGetLocalVecArray(iga,vecA,&localA,&arrayA);CHKERRQ(ierr);
  ierr = IGAGetLocalVecArray(iga,vecV,&localV,&arrayV);CHKERRQ(ierr);
  ierr = IGAGetLocalVecArray(iga,vecU,&localU,&arrayU);CHKERRQ(ierr);
  ierr = IGAGetLocalVecArray(iga,th->Vnet,&localVnet,&arrayVnet);CHKERRQ(ierr);   //Mesh velocity for the fluid



  ierr = PetscLogEventBegin(IGA_FormIFunction,iga,vecV,vecU,vecF);CHKERRQ(ierr);

  /* Element loop */
  ierr = IGABeginElement(iga,&element);CHKERRQ(ierr);
  while (IGANextElement(iga,element)) {


    if (iga->cruz == 1 || iga->cruz == 2 ) ierr = IGAElementFixValues0_JB(iga,element);CHKERRQ(ierr); //Added by me


    ierr = IGAElementGetWorkVal(element,&A);CHKERRQ(ierr);
    ierr = IGAElementGetWorkVal(element,&V);CHKERRQ(ierr);
    ierr = IGAElementGetWorkVal(element,&U);CHKERRQ(ierr);
    ierr = IGAElementGetWorkVec(element,&F);CHKERRQ(ierr);
    ierr = IGAElementGetWorkVec(element,&Vn);CHKERRQ(ierr);
    ierr = IGAElementGetValues(element,arrayA,A);CHKERRQ(ierr);
    ierr = IGAElementGetValues(element,arrayV,V);CHKERRQ(ierr);
    ierr = IGAElementGetValues(element,arrayU,U);CHKERRQ(ierr);
    ierr = IGAElementGetValues(element,arrayVnet,Vn);CHKERRQ(ierr);

    if (iga->cruz == 1 || iga->cruz == 2 ){
    ierr = IGAElementFixValues_JB(element,U,V,A);CHKERRQ(ierr);
    } else {
    ierr = IGAElementFixValues(element,U);CHKERRQ(ierr);}


    while (IGAElementNextUserIFunction(element,&IFunction,&ctx)) {

      /* Quadrature loop */
      ierr = IGAElementBeginPoint(element,&point);CHKERRQ(ierr);
      while (IGAElementNextPoint(element,point)) {
        ierr = IGAPointGetWorkVec(point,&R);CHKERRQ(ierr);
//        ierr = IFunction(point,dt,a,A,v,V,t,U,R,ctx);CHKERRQ(ierr);   //Original one
        ierr = IFunction(point,dt,a,A,v,V,t,U,R,Vn,ctx);CHKERRQ(ierr);   //Mine: Vn
        ierr = IGAPointAddVec(point,R,F);CHKERRQ(ierr);
      }
      ierr = IGAElementEndPoint(element,&point);CHKERRQ(ierr);
    }

   if (iga->cruz == 1 || iga->cruz == 2 ){
   ierr = IGAElementFixFunction_JB(element,F);CHKERRQ(ierr);
   }else{
   ierr = IGAElementFixFunction(element,F);CHKERRQ(ierr);
   }

    ierr = IGAElementAssembleVec(element,F,vecF);CHKERRQ(ierr);
  }
  ierr = IGAEndElement(iga,&element);CHKERRQ(ierr);

  ierr = PetscLogEventEnd(IGA_FormIFunction,iga,vecV,vecU,vecF);CHKERRQ(ierr);

  /* Restore local vectors V,U and arrays */
  ierr = IGARestoreLocalVecArray(iga,vecA,&localA,&arrayA);CHKERRQ(ierr);
  ierr = IGARestoreLocalVecArray(iga,vecV,&localV,&arrayV);CHKERRQ(ierr);
  ierr = IGARestoreLocalVecArray(iga,vecU,&localU,&arrayU);CHKERRQ(ierr);
  ierr = IGARestoreLocalVecArray(iga,th->Vnet,&localVnet,&arrayVnet);CHKERRQ(ierr);

  /* Assemble global vector F */
  ierr = VecAssemblyBegin(vecF);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(vecF);CHKERRQ(ierr);

//  if(iga->cruz == 2) ierr = VecView(vecF,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);


  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAComputeIJacobian2"
PetscErrorCode IGAComputeIJacobian2(IGA iga,PetscReal dt,
                                    PetscReal a,Vec vecA,
                                    PetscReal v,Vec vecV,
                                    PetscReal t,Vec vecU,
                                    Mat matJ, TS ts)  //TS ts
{
  Vec               localA;
  Vec               localV;
  Vec               localU;
  Vec				localVnet;
  const PetscScalar *arrayA;
  const PetscScalar *arrayV;
  const PetscScalar *arrayU;
  const PetscScalar *arrayVnet;
  IGAElement        element;
  IGAPoint          point;
  IGAUserIJacobian2 IJacobian;
  void              *ctx;
  PetscScalar       *A,*V,*U,*J,*K,*Vn;
  PetscErrorCode    ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidHeaderSpecific(vecA,VEC_CLASSID,4);
  PetscValidHeaderSpecific(vecV,VEC_CLASSID,6);
  PetscValidHeaderSpecific(vecU,VEC_CLASSID,8);
  PetscValidHeaderSpecific(matJ,MAT_CLASSID,9);
  IGACheckSetUp(iga,1);
  IGACheckUserOp(iga,1,IJacobian2);
  IJacobian = iga->userops->IJacobian2;
  ctx       = iga->userops->IJacCtx;
  TS_Alpha2     *th    = (TS_Alpha2*)ts->data;
//  AppCtx *user = (AppCtx *)ctx;
//  PetscPrintf(PETSC_COMM_WORLD,"IGAComputeIJacobian2 \n");



  /* Clear global matrix J*/
  ierr = MatZeroEntries(matJ);CHKERRQ(ierr);

  /* Get local vectors A,V,U and arrays */
  ierr = IGAGetLocalVecArray(iga,vecA,&localA,&arrayA);CHKERRQ(ierr);
  ierr = IGAGetLocalVecArray(iga,vecV,&localV,&arrayV);CHKERRQ(ierr);
  ierr = IGAGetLocalVecArray(iga,vecU,&localU,&arrayU);CHKERRQ(ierr);
  ierr = IGAGetLocalVecArray(iga,th->Vnet,&localVnet,&arrayVnet);CHKERRQ(ierr);

  ierr = PetscLogEventBegin(IGA_FormIJacobian,iga,vecV,vecU,matJ);CHKERRQ(ierr);

  /* Element Loop */
  ierr = IGABeginElement(iga,&element);CHKERRQ(ierr);
  while (IGANextElement(iga,element)) {


    if (iga->cruz == 1 || iga->cruz == 2 ) ierr = IGAElementFixValues0_JB(iga,element);CHKERRQ(ierr); //Added by me


    ierr = IGAElementGetWorkVal(element,&A);CHKERRQ(ierr);
    ierr = IGAElementGetWorkVal(element,&V);CHKERRQ(ierr);
    ierr = IGAElementGetWorkVal(element,&U);CHKERRQ(ierr);
    ierr = IGAElementGetWorkMat(element,&Vn);CHKERRQ(ierr);
    ierr = IGAElementGetWorkMat(element,&J);CHKERRQ(ierr);

    ierr = IGAElementGetValues(element,arrayA,A);CHKERRQ(ierr);
    ierr = IGAElementGetValues(element,arrayV,V);CHKERRQ(ierr);
    ierr = IGAElementGetValues(element,arrayU,U);CHKERRQ(ierr);
    ierr = IGAElementGetValues(element,arrayVnet,Vn);CHKERRQ(ierr);


    if (iga->cruz == 1 || iga->cruz == 2 ){
    ierr = IGAElementFixValues_JB(element,U,V,A);CHKERRQ(ierr);
    }else {    
    ierr = IGAElementFixValues(element,U);CHKERRQ(ierr);}


    while (IGAElementNextUserIJacobian(element,&IJacobian,&ctx)) {
      /* Quadrature loop */
      ierr = IGAElementBeginPoint(element,&point);CHKERRQ(ierr);
      while (IGAElementNextPoint(element,point)) {
        ierr = IGAPointGetWorkMat(point,&K);CHKERRQ(ierr);
//        ierr = IJacobian(point,dt,a,A,v,V,t,U,K,ctx);CHKERRQ(ierr);		//Original one
        ierr = IJacobian(point,dt,a,A,v,V,t,U,K,Vn,ctx);CHKERRQ(ierr);  //Mine   Vn
        ierr = IGAPointAddMat(point,K,J);CHKERRQ(ierr);
      }
      ierr = IGAElementEndPoint(element,&point);CHKERRQ(ierr);
    }


    if (iga->cruz == 1 || iga->cruz == 2 ){
    ierr = IGAElementFixJacobian_JB(element,J,ts);CHKERRQ(ierr);
    }else{
    ierr = IGAElementFixJacobian(element,J);CHKERRQ(ierr);
    }


    ierr = IGAElementAssembleMat(element,J,matJ);CHKERRQ(ierr);
  }
  ierr = IGAEndElement(iga,&element);CHKERRQ(ierr);


  ierr = PetscLogEventEnd(IGA_FormIJacobian,iga,vecV,vecU,matJ);CHKERRQ(ierr);

  /* Restore local vectors A,V,U and arrays */
  ierr = IGARestoreLocalVecArray(iga,vecA,&localA,&arrayA);CHKERRQ(ierr);
  ierr = IGARestoreLocalVecArray(iga,vecV,&localV,&arrayV);CHKERRQ(ierr);
  ierr = IGARestoreLocalVecArray(iga,vecU,&localU,&arrayU);CHKERRQ(ierr);
  ierr = IGARestoreLocalVecArray(iga,th->Vnet,&localVnet,&arrayVnet);CHKERRQ(ierr);

  /* Assemble global matrix J*/
  ierr = MatAssemblyBegin(matJ,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd  (matJ,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);


 // if(iga->cruz == 2) ierr = MatView(matJ,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);


  //#############################


  PetscInt nodesX  = iga->geom_lwidth[0], nodesY  = iga->geom_lwidth[1], nodesZ  = iga->geom_lwidth[2];
  PetscInt gnodesX = iga->geom_gwidth[0], gnodesY = iga->geom_gwidth[1];
  PetscInt dof = iga->dof;
  PetscInt dim = iga->dim;

if(iga->cruz == 2){

  PetscInt j,k,l,m;
  PetscInt npd;  ierr = IGAGetPropertyDim(iga,&npd);CHKERRQ(ierr);
  for(m=0;m<nodesZ;m++) {
  for(l=0;l<nodesY;l++) {
  for(k=0;k<nodesX;k++) {
	   if(iga->propertyA[(m*gnodesX*gnodesY + l*gnodesX+k)*npd] != 1.){
	    	 for (j=0; j<dof; j++) {
			 if(j == 0 || j == 1+dim || j == 2+dim){ ierr = MatSetValueLocal(matJ,(m*gnodesX*gnodesY + l*gnodesX+ k)*dof +j,(m*gnodesX*gnodesY + l*gnodesX+ k)*dof +j,(1.0/(th->scale_J)),INSERT_VALUES);CHKERRQ(ierr); }
	    	 }
	   }
//////////////////////////////////////// nsk_inter (dof = 5 but no energy equation)
//  	 for (j=0; j<dof; j++) {
//		 if(j == 2+dim){ ierr = MatSetValueLocal(matJ,(m*gnodesX*gnodesY + l*gnodesX+ k)*dof +j,(m*gnodesX*gnodesY + l*gnodesX+ k)*dof +j,(1.0/(th->scale_J)),INSERT_VALUES);CHKERRQ(ierr); }
//  	 }
////////////////////////////////////////

  }
  }
  }

//Just for the circumference and the square_ring

//  for(m=0;m<nodesZ;m++) {
//  for(l=0;l<nodesY;l++) {
//  for(k=0;k<nodesX;k++) {
//	   if(iga->propertyA[(m*gnodesX*gnodesY + l*gnodesX+k)*npd + 1] != 8.){
//	    	 for (j=0; j<dof; j++) {
//			 ierr = MatSetValueLocal(matJ,(m*gnodesX*gnodesY + l*gnodesX+ k)*dof +j,(m*gnodesX*gnodesY + l*gnodesX+ k)*dof +j,(1.0/(th->Alpha_f*th->Beta*dt*dt)),INSERT_VALUES);CHKERRQ(ierr);
//	    	 }
//	   }
//  }
//  }
//  }



 }

    ierr = MatAssemblyBegin(matJ,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd  (matJ,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

//   if(iga->cruz == 2) ierr = MatView(matJ,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);


    //#############################

  if (iga->cruz == 2){

   PetscInt k,l,m;
   PetscInt npd;  ierr = IGAGetPropertyDim(iga,&npd);CHKERRQ(ierr);
   PetscScalar *arrayV1;
   ierr  = VecGetArray(th->V1,&arrayV1);CHKERRQ(ierr);


//   for(i=0;i<Nnodes;i++) {
// 	  for(k=0;k<dof;k++) {
// 		    PetscPrintf(PETSC_COMM_WORLD,"i: %d, arrayV1: %e, arrayX1: %e \n",i,arrayV1[i*dof+k], arrayX1[i*dof+k]);
// 	  }
//   }


   for(m=0;m<nodesZ;m++) {
   for(l=0;l<nodesY;l++) {
   for(k=0;k<nodesX;k++) {
      	if(iga->propertyA[(m*gnodesX*gnodesY + l*gnodesX+k)*npd+1] == 8.)  arrayV1[(m*nodesX*nodesY + l*nodesX+ k)*dof + dim+2 ] = iga->T_wall; //323.5;
   }
   }
   }


//   for(i=0;i<Nnodes;i++) {
// 	  for(k=0;k<dof;k++) {
// 		    PetscPrintf(PETSC_COMM_WORLD,"i: %d, arrayV1 POST: %e \n",i, arrayV1[i*dof+k]);
// 	  }
//   }


   VecRestoreArray(th->V1,&arrayV1);CHKERRQ(ierr);
  }


  if (iga->cruz == 1){

  PetscInt j,k,l,m;
  PetscScalar *arrayX1, *arrayV1, *arrayA1;
  PetscScalar *arraydAs;
  ierr  = VecGetArray(th->X1,&arrayX1);CHKERRQ(ierr);
  ierr  = VecGetArray(th->V1,&arrayV1);CHKERRQ(ierr);
  ierr  = VecGetArray(th->A1,&arrayA1);CHKERRQ(ierr);
  ierr  = VecGetArray(th->dAs,&arraydAs);CHKERRQ(ierr);
  PetscInt npdG;  ierr = IGAGetPropertyDim(iga,&npdG);CHKERRQ(ierr);
 // ierr = VecGetSize(localU,&N);

//     for(i=0;i<15;i++) {
//   	  for(k=0;k<2;k++)   PetscPrintf(PETSC_COMM_WORLD,"i: %d, arrayX1: %e, arrayV1: %e, Jac arraydAs: %e \n",i, arrayX1[i*2+k],arrayV1[i*2+k],arraydAs[i*2+k]);
//     }

  	 for(m=0;m<nodesZ;m++) {
     for(l=0;l<nodesY;l++) {
     for(k=0;k<nodesX;k++) {
          for(j=0;j<dim;j++) {
        	  if(iga->propertyA[(m*gnodesX*gnodesY + l*gnodesX+k)*npdG+1] == 8.)  {
        		  arrayX1[(m*nodesX*nodesY + l*nodesX+ k)*dof +j] += th->Beta*dt*dt*arraydAs[(m*nodesX*nodesY + l*nodesX+ k)*dof +j];
        		  arrayV1[(m*nodesX*nodesY + l*nodesX+ k)*dof +j] += th->Gamma*dt*arraydAs[(m*nodesX*nodesY + l*nodesX+ k)*dof +j];
        		  arrayA1[(m*nodesX*nodesY + l*nodesX+ k)*dof +j] += arraydAs[(m*nodesX*nodesY + l*nodesX+ k)*dof +j];
        	  }
          }
     }
     }
  	 }

//     for(i=0;i<15;i++) {
//   	  for(k=0;k<2;k++) PetscPrintf(PETSC_COMM_WORLD,"i: %d, arrayX1 POST: %e, arrayV1 POST: %e, arraydAs POST: %e \n",i, arrayX1[i*2+k], arrayV1[i*2+k], arraydAs[i*2+k]);
//     }


  VecRestoreArray(th->A1,&arrayA1);CHKERRQ(ierr);
  VecRestoreArray(th->V1,&arrayV1);CHKERRQ(ierr);
  VecRestoreArray(th->X1,&arrayX1);CHKERRQ(ierr);
  VecRestoreArray(th->dAs,&arraydAs);CHKERRQ(ierr);

  }




  PetscFunctionReturn(0);
}

PETSC_EXTERN PetscErrorCode IGATSFormIFunction (TS,PetscReal,Vec,Vec,Vec,void*);
PETSC_EXTERN PetscErrorCode IGATSFormIJacobian (TS,PetscReal,Vec,Vec,PetscReal,Mat*,Mat*,MatStructure*,void*);
PETSC_EXTERN PetscErrorCode IGATSFormIFunction2(TS,PetscReal,Vec,Vec,Vec,Vec,void*);
PETSC_EXTERN PetscErrorCode IGATSFormIJacobian2(TS,PetscReal,Vec,Vec,Vec,PetscReal,PetscReal,Mat*,Mat*,MatStructure*,void*);

#undef  __FUNCT__
#define __FUNCT__ "IGATSFormIFunction2"
PetscErrorCode IGATSFormIFunction2(TS ts,PetscReal t,Vec U,Vec V,Vec A,Vec F,void *ctx)
{
  IGA            iga = (IGA)ctx;
  PetscReal      dt,a=0,v=0;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ts,TS_CLASSID,1);
  PetscValidHeaderSpecific(U,VEC_CLASSID,3);
  PetscValidHeaderSpecific(V,VEC_CLASSID,4);
  PetscValidHeaderSpecific(A,VEC_CLASSID,5);
  PetscValidHeaderSpecific(F,VEC_CLASSID,6);
  PetscValidHeaderSpecific(iga,IGA_CLASSID,7);
  ierr = TSGetTimeStep(ts,&dt);CHKERRQ(ierr);
  ierr = IGAComputeIFunction2(iga,dt,a,A,v,V,t,U,F,ts);CHKERRQ(ierr); //,ts
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGATSFormIJacobian2"
PetscErrorCode IGATSFormIJacobian2(TS ts,PetscReal t,Vec U,Vec V,Vec A,PetscReal shiftV,PetscReal shiftA,Mat *J,Mat *P,MatStructure *m,void *ctx)
{
  IGA            iga = (IGA)ctx;
  PetscReal      dt,a=shiftA,v=shiftV;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ts,TS_CLASSID,1);
  PetscValidHeaderSpecific(U,VEC_CLASSID,3);
  PetscValidHeaderSpecific(V,VEC_CLASSID,4);
  PetscValidHeaderSpecific(A,VEC_CLASSID,5);
  PetscValidPointer(J,8);
  PetscValidHeaderSpecific(*J,MAT_CLASSID,8);
  PetscValidPointer(P,9);
  PetscValidHeaderSpecific(*P,MAT_CLASSID,9);
  PetscValidPointer(m,10);
  PetscValidHeaderSpecific(iga,IGA_CLASSID,10);
  ierr = TSGetTimeStep(ts,&dt);CHKERRQ(ierr);
  ierr = IGAComputeIJacobian2(iga,dt,a,A,v,V,t,U,*P,ts);CHKERRQ(ierr); //,ts
  if (*J != *P) {
    ierr = MatAssemblyBegin(*J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(*J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  }
  *m = SAME_NONZERO_PATTERN;
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGACreateTS2"
/*@
   IGACreateTS2 - Creates a TS (time stepper) which uses the same
   communicators as the IGA.

   Logically collective on IGA

   Input Parameter:
.  iga - the IGA context

   Output Parameter:
.  ts - the TS

   Level: normal

.keywords: IGA, create, TS
@*/
PetscErrorCode IGACreateTS2(IGA iga, TS *ts)
{
  MPI_Comm       comm;
  Vec            U,V;
  Vec            F;
  Mat            J;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidPointer(ts,2);




  ierr = IGAGetComm(iga,&comm);CHKERRQ(ierr);
  ierr = TSCreate(comm,ts);CHKERRQ(ierr);
  ierr = TSSetType(*ts,TSALPHA2);CHKERRQ(ierr);
  ierr = PetscObjectCompose((PetscObject)*ts,"IGA",(PetscObject)iga);CHKERRQ(ierr);
  ierr = IGASetOptionsHandlerTS(*ts);CHKERRQ(ierr);

  ierr = IGACreateVec(iga,&U);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&V);CHKERRQ(ierr);
  ierr = TSSetSolution2(*ts,U,V);CHKERRQ(ierr);
  ierr = VecDestroy(&U);CHKERRQ(ierr);
  ierr = VecDestroy(&V);CHKERRQ(ierr);

  ierr = IGACreateVec(iga,&F);CHKERRQ(ierr);
  ierr = TSSetIFunction (*ts,F,IGATSFormIFunction ,iga);CHKERRQ(ierr);
  ierr = TSSetIFunction2(*ts,F,IGATSFormIFunction2,iga);CHKERRQ(ierr);
  ierr = VecDestroy(&F);CHKERRQ(ierr);

  ierr = IGACreateMat(iga,&J);CHKERRQ(ierr);
  ierr = TSSetIJacobian (*ts,J,J,IGATSFormIJacobian, iga);CHKERRQ(ierr);
  ierr = TSSetIJacobian2(*ts,J,J,IGATSFormIJacobian2,iga);CHKERRQ(ierr);
  ierr = MatDestroy(&J);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
