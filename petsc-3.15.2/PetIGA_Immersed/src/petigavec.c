#include "petiga.h"
#include <petsc-private/vecimpl.h>

#if PETSC_VERSION_LE(3,3,0)
EXTERN_C_BEGIN
#endif
PETSC_EXTERN PetscErrorCode VecView_MPI_DA(Vec,PetscViewer);
PETSC_EXTERN PetscErrorCode VecLoad_Default_DA(Vec,PetscViewer);
#if PETSC_VERSION_LE(3,3,0)
EXTERN_C_END
#endif

#if PETSC_VERSION_LE(3,3,0)
static PetscErrorCode VecSetLayout(Vec v,PetscLayout map)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(v,VEC_CLASSID,1);
  ierr = PetscLayoutReference(map,&v->map);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#endif

#if PETSC_VERSION_LE(3,3,0)
#undef  __FUNCT__
#define __FUNCT__ "VecGetDM"
static PetscErrorCode VecGetDM(Vec v,DM *dm)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(v,VEC_CLASSID,1);
  PetscValidPointer(dm,2);
  ierr = PetscObjectQuery((PetscObject)v,"DM",(PetscObject*)dm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#undef  __FUNCT__
#define __FUNCT__ "VecSetDM"
static PetscErrorCode VecSetDM(Vec v,DM dm)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(v,VEC_CLASSID,1);
  if (dm) PetscValidHeaderSpecific(dm,DM_CLASSID,2);
  ierr = PetscObjectCompose((PetscObject)v,"DM",(PetscObject)dm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#endif

#if PETSC_VERSION_(3,4,0)
#define VecSetDM(v,dm) PetscObjectCompose((PetscObject)v,"__PETSc_dm",(PetscObject)dm)
#endif


#undef  __FUNCT__
#define __FUNCT__ "VecDuplicate_IGA"
static PetscErrorCode VecDuplicate_IGA(Vec g,Vec* gg)
{
  IGA            iga;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = PetscObjectQuery((PetscObject)g,"IGA",(PetscObject*)&iga);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,gg);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "VecView_IGA"
static PetscErrorCode VecView_IGA(Vec v,PetscViewer viewer)
{
  IGA            iga;
  DM             save,dm;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = VecGetDM(v,&save);CHKERRQ(ierr);
  if (save) {ierr = PetscObjectReference((PetscObject)save);CHKERRQ(ierr);}
  ierr = PetscObjectQuery((PetscObject)v,"IGA",(PetscObject*)&iga);CHKERRQ(ierr);
  ierr = IGAGetNodeDM(iga,&dm);CHKERRQ(ierr);
  ierr = VecSetDM(v,dm);CHKERRQ(ierr);
  ierr = VecView_MPI_DA(v,viewer);CHKERRQ(ierr);
  ierr = VecSetDM(v,save);CHKERRQ(ierr);
  ierr = DMDestroy(&save);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#undef  __FUNCT__
#define __FUNCT__ "VecLoad_IGA"
static PetscErrorCode VecLoad_IGA(Vec v,PetscViewer viewer)
{
  IGA            iga;
  DM             save,dm;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = VecGetDM(v,&save);CHKERRQ(ierr);
  if (save) {ierr = PetscObjectReference((PetscObject)save);CHKERRQ(ierr);}
  ierr = PetscObjectQuery((PetscObject)v,"IGA",(PetscObject*)&iga);CHKERRQ(ierr);
  ierr = IGAGetNodeDM(iga,&dm);CHKERRQ(ierr);
  ierr = VecSetDM(v,dm);CHKERRQ(ierr);
  ierr = VecLoad_Default_DA(v,viewer);CHKERRQ(ierr);
  ierr = VecSetDM(v,save);CHKERRQ(ierr);
  ierr = DMDestroy(&save);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PETSC_STATIC_INLINE
PetscInt Product(const PetscInt a[3]) { return a[0]*a[1]*a[2]; }

#undef  __FUNCT__
#define __FUNCT__ "IGACreateVec"
/*@
   IGACreateVec - Creates a vector with the correct parallel layout
   required for computing a vector using the discretization
   information provided in the IGA.

   Collective on IGA

   Input Parameter:
.  iga - the IGA context

   Output Parameter:
.  vec - the vector

   Level: normal

.keywords: IGA, create, vector
@*/
PetscErrorCode IGACreateVec(IGA iga, Vec *vec)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidPointer(vec,2);
  IGACheckSetUp(iga,1);
  ierr = VecCreate(((PetscObject)iga)->comm,vec);CHKERRQ(ierr);
  ierr = VecSetLayout(*vec,iga->map);CHKERRQ(ierr);
  ierr = VecSetType(*vec,iga->vectype);CHKERRQ(ierr);
  ierr = VecSetFromOptions(*vec);CHKERRQ(ierr);
  ierr = PetscObjectCompose((PetscObject)*vec,"IGA",(PetscObject)iga);CHKERRQ(ierr);
  ierr = VecSetOperation(*vec,VECOP_DUPLICATE,(void(*)(void))VecDuplicate_IGA);CHKERRQ(ierr);
  ierr = VecSetOperation(*vec,VECOP_VIEW,(void(*)(void))VecView_IGA);CHKERRQ(ierr);
  ierr = VecSetOperation(*vec,VECOP_LOAD,(void(*)(void))VecLoad_IGA);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGACreateLocalVec"
PetscErrorCode IGACreateLocalVec(IGA iga, Vec *vec)
{
  PetscInt       bs,n;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidPointer(vec,2);
  IGACheckSetUp(iga,1);
  /* */
  bs = iga->dof;
  n  = Product(iga->node_gwidth);
  ierr = VecCreate(PETSC_COMM_SELF,vec);CHKERRQ(ierr);
  ierr = VecSetSizes(*vec,n*bs,n*bs);CHKERRQ(ierr);
  ierr = VecSetBlockSize(*vec,bs);CHKERRQ(ierr);
  ierr = VecSetType(*vec,iga->vectype);CHKERRQ(ierr);
  ierr = VecSetFromOptions(*vec);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAGetLocalVec"
PetscErrorCode IGAGetLocalVec(IGA iga,Vec *lvec)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidPointer(lvec,2);
  IGACheckSetUp(iga,1);
  if (iga->nwork > 0) {
    *lvec = iga->vwork[--iga->nwork];
    iga->vwork[iga->nwork] = NULL;
  } else {
    ierr = IGACreateLocalVec(iga,lvec);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGARestoreLocalVec"
PetscErrorCode IGARestoreLocalVec(IGA iga,Vec *lvec)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidPointer(lvec,2);
  PetscValidHeaderSpecific(*lvec,VEC_CLASSID,2);
  IGACheckSetUp(iga,1);
  if (iga->nwork < (PetscInt)(sizeof(iga->vwork)/sizeof(Vec))) {
    iga->vwork[iga->nwork++] = *lvec; *lvec = 0;
  } else {
    ierr = VecDestroy(lvec);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAGlobalToLocalBegin"
PetscErrorCode IGAGlobalToLocalBegin(IGA iga,Vec gvec,Vec lvec,InsertMode addv)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidHeaderSpecific(gvec,VEC_CLASSID,2);
  PetscValidHeaderSpecific(lvec,VEC_CLASSID,3);
  IGACheckSetUp(iga,1);
  ierr = VecScatterBegin(iga->g2l,gvec,lvec,addv,SCATTER_FORWARD);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAGlobalToLocalEnd"
PetscErrorCode IGAGlobalToLocalEnd(IGA iga,Vec gvec,Vec lvec,InsertMode addv)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidHeaderSpecific(gvec,VEC_CLASSID,2);
  PetscValidHeaderSpecific(lvec,VEC_CLASSID,3);
  IGACheckSetUp(iga,1);
  ierr = VecScatterEnd(iga->g2l,gvec,lvec,addv,SCATTER_FORWARD);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAGlobalToLocal"
PetscErrorCode IGAGlobalToLocal(IGA iga,Vec gvec,Vec lvec,InsertMode addv)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = IGAGlobalToLocalBegin(iga,gvec,lvec,addv);CHKERRQ(ierr);
  ierr = IGAGlobalToLocalEnd  (iga,gvec,lvec,addv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGALocalToGlobalBegin"
PetscErrorCode IGALocalToGlobalBegin(IGA iga,Vec lvec,Vec gvec,InsertMode addv)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidHeaderSpecific(lvec,VEC_CLASSID,2);
  PetscValidHeaderSpecific(gvec,VEC_CLASSID,3);
  IGACheckSetUp(iga,1);
  if (addv == ADD_VALUES) {
    ierr = VecScatterBegin(iga->g2l,lvec,gvec,ADD_VALUES,SCATTER_REVERSE);CHKERRQ(ierr);
  } else if (addv == INSERT_VALUES) {
    ierr = VecScatterBegin(iga->l2g,lvec,gvec,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  } else SETERRQ(((PetscObject)iga)->comm,PETSC_ERR_SUP,"Not yet implemented");
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGALocalToGlobalEnd"
PetscErrorCode IGALocalToGlobalEnd(IGA iga,Vec lvec,Vec gvec,InsertMode addv)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidHeaderSpecific(lvec,VEC_CLASSID,2);
  PetscValidHeaderSpecific(gvec,VEC_CLASSID,3);
  IGACheckSetUp(iga,1);
  if (addv == ADD_VALUES) {
    ierr = VecScatterEnd(iga->g2l,lvec,gvec,ADD_VALUES,SCATTER_REVERSE);CHKERRQ(ierr);
  } else if (addv == INSERT_VALUES) {
    ierr = VecScatterEnd(iga->l2g,lvec,gvec,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  } else SETERRQ(((PetscObject)iga)->comm,PETSC_ERR_SUP,"Not yet implemented");
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGALocalToGlobal"
PetscErrorCode IGALocalToGlobal(IGA iga,Vec lvec,Vec gvec,InsertMode addv)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = IGALocalToGlobalBegin(iga,lvec,gvec,addv);CHKERRQ(ierr);
  ierr = IGALocalToGlobalEnd  (iga,lvec,gvec,addv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAGetLocalVecArray"
PetscErrorCode IGAGetLocalVecArray(IGA iga,Vec gvec,Vec *lvec,const PetscScalar *array[])
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidHeaderSpecific(gvec,VEC_CLASSID,2);
  PetscValidPointer(lvec,3);
  PetscValidPointer(array,4);
  IGACheckSetUp(iga,1);
  ierr = IGAGetLocalVec(iga,lvec);CHKERRQ(ierr);
  ierr = IGAGlobalToLocal(iga,gvec,*lvec,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecGetArrayRead(*lvec,array);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAGetLocalVecArrayModified"
PetscErrorCode IGAGetLocalVecArrayModified(IGA iga,Vec gvec,Vec *lvec,PetscScalar *array[])
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
//  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
//  PetscValidHeaderSpecific(gvec,VEC_CLASSID,2);
//  PetscValidPointer(lvec,3);
//  PetscValidPointer(array,4);
//  IGACheckSetUp(iga,1);
  ierr = IGAGetLocalVec(iga,lvec);CHKERRQ(ierr);
  ierr = IGAGlobalToLocal(iga,gvec,*lvec,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecGetArray(*lvec,array);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGARestoreLocalVecArray"
PetscErrorCode IGARestoreLocalVecArray(IGA iga,Vec gvec,Vec *lvec,const PetscScalar *array[])
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidHeaderSpecific(gvec,VEC_CLASSID,2);
  PetscValidPointer(lvec,3);
  PetscValidHeaderSpecific(*lvec,VEC_CLASSID,3);
  PetscValidPointer(array,4);
  IGACheckSetUp(iga,1);
  ierr = VecRestoreArrayRead(*lvec,array);CHKERRQ(ierr);
  ierr = IGARestoreLocalVec(iga,lvec);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGARestoreLocalVecArrayModified"
PetscErrorCode IGARestoreLocalVecArrayModified(IGA iga,Vec gvec,Vec *lvec,PetscScalar *array[])
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
//  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
//  PetscValidHeaderSpecific(gvec,VEC_CLASSID,2);
//  PetscValidPointer(lvec,3);
//  PetscValidHeaderSpecific(*lvec,VEC_CLASSID,3);
//  PetscValidPointer(array,4);
//  IGACheckSetUp(iga,1);
  ierr = VecRestoreArray(*lvec,array);CHKERRQ(ierr);
  ierr = IGARestoreLocalVec(iga,lvec);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#undef  __FUNCT__
#define __FUNCT__ "IGAGetNaturalVec"
PetscErrorCode IGAGetNaturalVec(IGA iga,Vec *nvec)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidPointer(nvec,2);
  IGACheckSetUp(iga,1);
  *nvec = iga->natural;
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGANaturalToGlobal"
PetscErrorCode IGANaturalToGlobal(IGA iga,Vec nvec,Vec gvec)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidHeaderSpecific(nvec,VEC_CLASSID,2);
  PetscValidHeaderSpecific(gvec,VEC_CLASSID,3);
  IGACheckSetUp(iga,1);
  ierr = VecScatterBegin(iga->n2g,nvec,gvec,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  ierr = VecScatterEnd  (iga->n2g,nvec,gvec,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAGlobalToNatural"
PetscErrorCode IGAGlobalToNatural(IGA iga,Vec gvec,Vec nvec)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidHeaderSpecific(gvec,VEC_CLASSID,2);
  PetscValidHeaderSpecific(nvec,VEC_CLASSID,3);
  IGACheckSetUp(iga,1);
  ierr = VecScatterBegin(iga->g2n,gvec,nvec,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  ierr = VecScatterEnd  (iga->g2n,gvec,nvec,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
