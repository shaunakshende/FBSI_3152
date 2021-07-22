#include "petiga.h"

#undef  __FUNCT__
#define __FUNCT__ "IGAInterfaceCreate"
PetscErrorCode IGAInterfaceCreate(IGAInterface *interface)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidPointer(interface,1);
  ierr = PetscNew(struct _n_IGAInterface,interface);CHKERRQ(ierr);
  (*interface)->refct = 1;
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAInterfaceDestroy"
PetscErrorCode IGAInterfaceDestroy(IGAInterface *_interface)
{
  IGAInterface    interface;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidPointer(_interface,1);
  interface = *_interface; *_interface = 0;
  if (!interface) PetscFunctionReturn(0);
  if (--interface->refct > 0) PetscFunctionReturn(0);
  ierr = IGAInterfaceReset(interface);CHKERRQ(ierr);
  ierr = PetscFree(interface);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAInterfaceReset"
PetscErrorCode IGAInterfaceReset(IGAInterface interface)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if (!interface) PetscFunctionReturn(0);
  PetscValidPointer(interface,1);
  interface->dof = 0;
  interface->count = 0;
  ierr = PetscFree(interface->field);CHKERRQ(ierr);
  ierr = PetscFree(interface->value);CHKERRQ(ierr);
  interface->nload = 0;
  ierr = PetscFree(interface->iload);CHKERRQ(ierr);
  ierr = PetscFree(interface->vload);CHKERRQ(ierr);
  ierr = PetscFree(interface->userops);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAInterfaceReference"
PetscErrorCode IGAInterfaceReference(IGAInterface interface)
{
  PetscFunctionBegin;
  PetscValidPointer(interface,1);
  interface->refct++;
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAInterfaceInit"
PetscErrorCode IGAInterfaceInit(IGAInterface interface,PetscInt dof)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidPointer(interface,1);
  ierr = IGAInterfaceReset(interface);CHKERRQ(ierr);
  interface->dof = dof;
  interface->count = 0;
  ierr = PetscMalloc1(dof,PetscInt,   &interface->field);CHKERRQ(ierr);
  ierr = PetscMalloc1(dof,PetscScalar,&interface->value);CHKERRQ(ierr);
  interface->nload = 0;
  ierr = PetscMalloc1(dof,PetscInt,   &interface->iload);CHKERRQ(ierr);
  ierr = PetscMalloc1(dof,PetscScalar,&interface->vload);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAInterfaceClear"
PetscErrorCode IGAInterfaceClear(IGAInterface interface)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidPointer(interface,1);
  interface->count = 0;
  interface->nload = 0;
  ierr = PetscFree(interface->userops);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAInterfaceSetValue"
/*@
   IGAInterfaceSetValue - Used to set a constant Dirichlet interface
   condition on the given interface.

   Logically Collective on IGAInterface

   Input Parameters:
+  interface - the IGAAxis context
.  field - the index of the field on which to enforce the condition
-  value - the value to set

   Level: normal

.keywords: IGA, interface, Dirichlet
@*/
PetscErrorCode IGAInterfaceSetValue(IGAInterface interface,PetscInt field,PetscScalar value)
{
  PetscInt dof;
  PetscFunctionBegin;
  PetscValidPointer(interface,1);
  dof = interface->dof;
  PetscPrintf(PETSC_COMM_WORLD,"Interface set Value (Dirichlet). Init!\n");
  if (field <  0)   SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Field %D must be nonnegative",field);
  if (field >= dof) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Field %D, but dof %D",field,dof);
  { /**/
    PetscInt pos;
    for (pos=0; pos<interface->count; pos++)
      if (interface->field[pos] == field) break;
    if (pos==interface->count) interface->count++;
    interface->field[pos] = field;
    interface->value[pos] = value;
  }
  PetscFunctionReturn(0);
  PetscPrintf(PETSC_COMM_WORLD,"Interface set Value (Dirichlet). End!\n");
}

#undef  __FUNCT__
#define __FUNCT__ "IGAInterfaceSetLoad"
/*@
   IGAInterfaceSetLoad - Used to set a constant Neumann interface
   condition on the given interface.

   Logically Collective on IGAInterface

   Input Parameters:
+  interface - the IGAAxis context
.  field - the index of the field on which to enforce the condition
-  value - the value to set

   Level: normal

.keywords: IGA, interface, Neumann
@*/
PetscErrorCode IGAInterfaceSetLoad(IGAInterface interface,PetscInt field,PetscScalar value)
{
  PetscInt dof;
  PetscFunctionBegin;
  PetscValidPointer(interface,1);
  dof = interface->dof;
  if (field <  0)   SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Field %D must be nonnegative",field);
  if (field >= dof) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Field %D, but dof %D",field,dof);
  { /**/
    PetscInt pos;
    for (pos=0; pos<interface->nload; pos++)
      if (interface->iload[pos] == field) break;
    if (pos==interface->nload) interface->nload++;
    interface->iload[pos] = field;
    interface->vload[pos] = value;
  }
  PetscFunctionReturn(0);
}

#define IGAInterfaceEnsureUserOps(interface) do { PetscErrorCode _ierr;                                     \
    if (!(interface)->userops) {_ierr = PetscNew(struct _IGAUserOps,&(interface)->userops);CHKERRQ(_ierr);} \
  } while (0)

#undef  __FUNCT__
#define __FUNCT__ "IGAInterfaceSetUserVector"
PetscErrorCode IGAInterfaceSetUserVector(IGAInterface interface,IGAUserVector Vector,void *VecCtx)
{
  PetscFunctionBegin;
  PetscValidPointer(interface,1);
  IGAInterfaceEnsureUserOps(interface);
  if (Vector) interface->userops->Vector = Vector;
  if (VecCtx) interface->userops->VecCtx = VecCtx;
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAInterfaceSetUserMatrix"
PetscErrorCode IGAInterfaceSetUserMatrix(IGAInterface interface,IGAUserMatrix Matrix,void *MatCtx)
{
  PetscFunctionBegin;
  PetscValidPointer(interface,1);
  IGAInterfaceEnsureUserOps(interface);
  if (Matrix) interface->userops->Matrix = Matrix;
  if (MatCtx) interface->userops->MatCtx = MatCtx;
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAInterfaceSetUserSystem"
/*@
   IGAInterfaceSetUserSystem - Set the user callback to form the matrix
   and vector which represents the discretized a(w,u) = L(w)
   integrated along the given interface.

   Logically collective on IGAInterface

   Input Parameters:
+  interface - the IGAInterface context
.  System - the function which evaluates a(w,u) and L(w)
-  ctx - user-defined context for evaluation routine (may be NULL)

   Details of System:
$  PetscErrorCode System(IGAPoint p,PetscScalar *K,PetscScalar *F,void *ctx);

+  p - point at which to evaluate a(w,u)=L(w)
.  K - contribution to a(w,u)
.  F - contribution to L(w)
-  ctx - user-defined context for evaluation routine

   Level: normal

.keywords: IGAInterface, setup linear system, matrix assembly, vector assembly
@*/
PetscErrorCode IGAInterfaceSetUserSystem(IGAInterface interface,IGAUserSystem System,void *SysCtx)
{
  PetscFunctionBegin;
  PetscValidPointer(interface,1);
  IGAInterfaceEnsureUserOps(interface);
  if (System) interface->userops->System = System;
  if (SysCtx) interface->userops->SysCtx = SysCtx;
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAInterfaceSetUserFunction"
PetscErrorCode IGAInterfaceSetUserFunction(IGAInterface interface,IGAUserFunction Function,void *FunCtx)
{
  PetscFunctionBegin;
  PetscValidPointer(interface->userops,0);
  IGAInterfaceEnsureUserOps(interface);
  if (Function) interface->userops->Function = Function;
  if (FunCtx)   interface->userops->FunCtx   = FunCtx;
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAInterfaceSetUserJacobian"
PetscErrorCode IGAInterfaceSetUserJacobian(IGAInterface interface,IGAUserJacobian Jacobian,void *JacCtx)
{
  PetscFunctionBegin;
  PetscValidPointer(interface,1);
  IGAInterfaceEnsureUserOps(interface);
  if (Jacobian) interface->userops->Jacobian = Jacobian;
  if (JacCtx)   interface->userops->JacCtx   = JacCtx;
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAInterfaceSetUserIFunction"
PetscErrorCode IGAInterfaceSetUserIFunction(IGAInterface interface,IGAUserIFunction IFunction,void *IFunCtx)
{
  PetscFunctionBegin;
  PetscValidPointer(interface,1);
  IGAInterfaceEnsureUserOps(interface);
  if (IFunction) interface->userops->IFunction = IFunction;
  if (IFunCtx)   interface->userops->IFunCtx   = IFunCtx;
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAInterfaceSetUserIJacobian"
PetscErrorCode IGAInterfaceSetUserIJacobian(IGAInterface interface,IGAUserIJacobian IJacobian,void *IJacCtx)
{
  PetscFunctionBegin;
  PetscValidPointer(interface,1);
  IGAInterfaceEnsureUserOps(interface);
  if (IJacobian) interface->userops->IJacobian = IJacobian;
  if (IJacCtx)   interface->userops->IJacCtx   = IJacCtx;
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAInterfaceSetUserIFunction2"
PetscErrorCode IGAInterfaceSetUserIFunction2(IGAInterface interface,IGAUserIFunction2 IFunction,void *IFunCtx)
{
  PetscFunctionBegin;
  PetscValidPointer(interface,1);
  IGAInterfaceEnsureUserOps(interface);
  if (IFunction) interface->userops->IFunction2 = IFunction;
  if (IFunCtx)   interface->userops->IFunCtx    = IFunCtx;
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAInterfaceSetUserIJacobian2"
PetscErrorCode IGAInterfaceSetUserIJacobian2(IGAInterface interface,IGAUserIJacobian2 IJacobian,void *IJacCtx)
{
  PetscFunctionBegin;
  PetscValidPointer(interface,1);
  IGAInterfaceEnsureUserOps(interface);
  if (IJacobian) interface->userops->IJacobian2 = IJacobian;
  if (IJacCtx)   interface->userops->IJacCtx    = IJacCtx;
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAInterfaceSetUserIEFunction"
PetscErrorCode IGAInterfaceSetUserIEFunction(IGAInterface interface,IGAUserIEFunction IEFunction,void *IEFunCtx)
{
  PetscFunctionBegin;
  PetscValidPointer(interface,1);
  IGAInterfaceEnsureUserOps(interface);
  if (IEFunction) interface->userops->IEFunction = IEFunction;
  if (IEFunCtx)   interface->userops->IEFunCtx   = IEFunCtx;
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAInterfaceSetUserIEJacobian"
PetscErrorCode IGAInterfaceSetUserIEJacobian(IGAInterface interface,IGAUserIEJacobian IEJacobian,void *IEJacCtx)
{
  PetscFunctionBegin;
  PetscValidPointer(interface,1);
  IGAInterfaceEnsureUserOps(interface);
  if (IEJacobian) interface->userops->IEJacobian = IEJacobian;
  if (IEJacCtx)   interface->userops->IEJacCtx   = IEJacCtx;
  PetscFunctionReturn(0);
}
