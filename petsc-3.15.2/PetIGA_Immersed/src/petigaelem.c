#include "petiga.h"
#include <petsc-private/tsimpl.h>
#include "jb.h"

#undef  __FUNCT__
#define __FUNCT__ "IGAElementCreate"
PetscErrorCode IGAElementCreate(IGAElement *_element)
{
  IGAElement element;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidPointer(_element,1);
  ierr = PetscNew(struct _n_IGAElement,_element);CHKERRQ(ierr);
  element = *_element;
  element->refct =  1;
  element->index = -1;
  ierr = IGAPointCreate(&element->iterator);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAElementDestroy"
PetscErrorCode IGAElementDestroy(IGAElement *_element)
{
  IGAElement     element;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidPointer(_element,1);
  element = *_element; *_element = 0;
  if (!element) PetscFunctionReturn(0);
  if (--element->refct > 0) PetscFunctionReturn(0);
  ierr = IGAPointDestroy(&element->iterator);CHKERRQ(ierr);
  ierr = IGAElementReset(element);CHKERRQ(ierr);
  ierr = PetscFree(element);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAElementReference"
PetscErrorCode IGAElementReference(IGAElement element)
{
  PetscFunctionBegin;
  PetscValidPointer(element,1);
  element->refct++;
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAElementFreeWork"
static
PetscErrorCode IGAElementFreeWork(IGAElement element)
{
  size_t i;
  size_t MAX_WORK_VAL;
  size_t MAX_WORK_VEC;
  size_t MAX_WORK_MAT;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidPointer(element,1);
  MAX_WORK_VAL = sizeof(element->wval)/sizeof(PetscScalar*);
  for (i=0; i<MAX_WORK_VAL; i++)
    {ierr = PetscFree(element->wval[i]);CHKERRQ(ierr);}
  element->nval = 0;
  MAX_WORK_VEC = sizeof(element->wvec)/sizeof(PetscScalar*);
  for (i=0; i<MAX_WORK_VEC; i++)
    {ierr = PetscFree(element->wvec[i]);CHKERRQ(ierr);}
  element->nvec = 0;
  MAX_WORK_MAT = sizeof(element->wmat)/sizeof(PetscScalar*);
  for (i=0; i<MAX_WORK_MAT; i++)
    {ierr = PetscFree(element->wmat[i]);CHKERRQ(ierr);}
  element->nmat = 0;
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAElementFreeFix"
static
PetscErrorCode IGAElementFreeFix(IGAElement element)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidPointer(element,1);
  element->nfix = 0;
  ierr = PetscFree(element->ifix);CHKERRQ(ierr);
  ierr = PetscFree(element->vfix);CHKERRQ(ierr);
  ierr = PetscFree(element->ufix);CHKERRQ(ierr);
  element->nflux = 0;
  ierr = PetscFree(element->iflux);CHKERRQ(ierr);
  ierr = PetscFree(element->vflux);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAElementReset"
PetscErrorCode IGAElementReset(IGAElement element)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if (!element) PetscFunctionReturn(0);
  PetscValidPointer(element,1);
  element->count =  0;
  element->index = -1;

  if (element->rowmap != element->mapping)
    {ierr = PetscFree(element->rowmap);CHKERRQ(ierr);}
  if (element->colmap != element->mapping)
    {ierr = PetscFree(element->colmap);CHKERRQ(ierr);}
  element->rowmap = element->colmap = 0;
  ierr = PetscFree(element->mapping);CHKERRQ(ierr);

  ierr = PetscFree(element->geometryX);CHKERRQ(ierr);
  ierr = PetscFree(element->geometryINI);CHKERRQ(ierr);
  ierr = PetscFree(element->rationalW);CHKERRQ(ierr);
  ierr = PetscFree(element->propertyA);CHKERRQ(ierr);

  ierr = PetscFree(element->weight);CHKERRQ(ierr);
  ierr = PetscFree(element->detJac);CHKERRQ(ierr);

  ierr = PetscFree(element->point);CHKERRQ(ierr);
  ierr = PetscFree(element->basis[0]);CHKERRQ(ierr);
  ierr = PetscFree(element->basis[1]);CHKERRQ(ierr);
  ierr = PetscFree(element->basis[2]);CHKERRQ(ierr);
  ierr = PetscFree(element->basis[3]);CHKERRQ(ierr);

  ierr = PetscFree(element->detX);CHKERRQ(ierr);
  ierr = PetscFree(element->gradX[0]);CHKERRQ(ierr);
  ierr = PetscFree(element->gradX[1]);CHKERRQ(ierr);
  ierr = PetscFree(element->shape[0]);CHKERRQ(ierr);
  ierr = PetscFree(element->shape[1]);CHKERRQ(ierr);
  ierr = PetscFree(element->shape[2]);CHKERRQ(ierr);
  ierr = PetscFree(element->shape[3]);CHKERRQ(ierr);

  ierr = PetscFree(element->normal);CHKERRQ(ierr);

  ierr = IGAElementFreeFix(element);CHKERRQ(ierr);
  ierr = IGAElementFreeWork(element);CHKERRQ(ierr);
  ierr = IGAPointReset(element->iterator);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAElementInit"
PetscErrorCode IGAElementInit(IGAElement element,IGA iga)
{
  PetscInt *start;
  PetscInt *width;
  PetscInt *sizes;
  IGABasis *BD;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidPointer(element,1);
  PetscValidHeaderSpecific(iga,IGA_CLASSID,2);
  IGACheckSetUp(iga,1);
  ierr = IGAElementReset(element);CHKERRQ(ierr);
  element->parent = iga;
  element->collocation = iga->collocation;

  element->dof = iga->dof;
  element->dim = iga->dim;
  element->nsd = iga->geometry ? iga->geometry : iga->dim;
  element->npd = iga->property ? iga->property : 0;

  start = iga->elem_start;
  width = iga->elem_width;
  sizes = iga->elem_sizes;
  BD    = iga->basis;

  { /* */
    PetscInt i,dim = element->dim;
    PetscInt nel=1,nen=1,nqp=1;
    for (i=0; i<3; i++) {
      element->ID[i] = 0;
      element->BD[i] = BD[i];
    }
    for (i=0; i<dim; i++) {
      element->start[i] = start[i];
      element->width[i] = width[i];
      element->sizes[i] = sizes[i];
      nel *= width[i];
      nen *= BD[i]->nen;
      nqp *= BD[i]->nqp;
    }
    for (i=dim; i<3; i++) {
      element->start[i] = 0;
      element->width[i] = 1;
      element->sizes[i] = 1;
    }
    element->index = -1;
    element->count = nel;
    element->nen   = nen;
    element->nqp   = nqp;
  }
  { /**/
    PetscInt nen = element->nen;
    ierr = PetscMalloc1(nen,PetscInt,&element->mapping);CHKERRQ(ierr);
    if (!element->collocation) {
      element->neq = nen;
      element->rowmap = element->mapping;
    } else {
      element->neq = 1;
      ierr = PetscMalloc1(1,PetscInt,&element->rowmap);CHKERRQ(ierr);
    }
    element->colmap = element->mapping;
  }
  { /* */
    PetscInt dim = element->dim;
    PetscInt nsd = element->nsd;
    PetscInt npd = element->npd;
    PetscInt nen = element->nen;
    PetscInt nqp = element->nqp;

    /* */

    ierr = PetscMalloc1(nen*nsd,PetscReal,&element->geometryX);CHKERRQ(ierr);
    ierr = PetscMalloc1(nen*nsd,PetscReal,&element->geometryINI);CHKERRQ(ierr);
    ierr = PetscMalloc1(nen,PetscReal,&element->rationalW);CHKERRQ(ierr);
    ierr = PetscMalloc1(nen*npd,PetscScalar,&element->propertyA);CHKERRQ(ierr);

    ierr = PetscMalloc1(nqp,PetscReal,&element->weight);CHKERRQ(ierr);
    ierr = PetscMalloc1(nqp,PetscReal,&element->detJac);CHKERRQ(ierr);

    ierr = PetscMalloc1(nqp*dim,PetscReal,&element->point);CHKERRQ(ierr);
    ierr = PetscMalloc1(nqp*nen,PetscReal,&element->basis[0]);CHKERRQ(ierr);
    ierr = PetscMalloc1(nqp*nen*dim,PetscReal,&element->basis[1]);CHKERRQ(ierr);
    ierr = PetscMalloc1(nqp*nen*dim*dim,PetscReal,&element->basis[2]);CHKERRQ(ierr);
    ierr = PetscMalloc1(nqp*nen*dim*dim*dim,PetscReal,&element->basis[3]);CHKERRQ(ierr);

    ierr = PetscMalloc1(nqp,PetscReal,&element->detX);CHKERRQ(ierr);
    ierr = PetscMalloc1(nqp*dim*dim,PetscReal,&element->gradX[0]);CHKERRQ(ierr);
    ierr = PetscMalloc1(nqp*dim*dim,PetscReal,&element->gradX[1]);CHKERRQ(ierr);
    ierr = PetscMalloc1(nqp*nen,PetscReal,&element->shape[0]);CHKERRQ(ierr);
    ierr = PetscMalloc1(nqp*nen*dim,PetscReal,&element->shape[1]);CHKERRQ(ierr);
    ierr = PetscMalloc1(nqp*nen*dim*dim,PetscReal,&element->shape[2]);CHKERRQ(ierr);
    ierr = PetscMalloc1(nqp*nen*dim*dim*dim,PetscReal,&element->shape[3]);CHKERRQ(ierr);

    ierr = PetscMalloc1(nqp*dim,PetscReal,&element->normal);CHKERRQ(ierr);

    /* */

    ierr = PetscMemzero(element->weight,  sizeof(PetscReal)*nqp);CHKERRQ(ierr);
    ierr = PetscMemzero(element->detJac,  sizeof(PetscReal)*nqp);CHKERRQ(ierr);

    ierr = PetscMemzero(element->point,   sizeof(PetscReal)*nqp*dim);CHKERRQ(ierr);
    ierr = PetscMemzero(element->basis[0],sizeof(PetscReal)*nqp*nen);CHKERRQ(ierr);
    ierr = PetscMemzero(element->basis[1],sizeof(PetscReal)*nqp*nen*dim);CHKERRQ(ierr);
    ierr = PetscMemzero(element->basis[2],sizeof(PetscReal)*nqp*nen*dim*dim);CHKERRQ(ierr);
    ierr = PetscMemzero(element->basis[3],sizeof(PetscReal)*nqp*nen*dim*dim*dim);CHKERRQ(ierr);

    ierr = PetscMemzero(element->detX,    sizeof(PetscReal)*nqp);CHKERRQ(ierr);
    ierr = PetscMemzero(element->gradX[0],sizeof(PetscReal)*nqp*dim*dim);CHKERRQ(ierr);
    ierr = PetscMemzero(element->gradX[1],sizeof(PetscReal)*nqp*dim*dim);CHKERRQ(ierr);
    ierr = PetscMemzero(element->shape[0],sizeof(PetscReal)*nqp*nen);CHKERRQ(ierr);
    ierr = PetscMemzero(element->shape[1],sizeof(PetscReal)*nqp*nen*dim);CHKERRQ(ierr);
    ierr = PetscMemzero(element->shape[2],sizeof(PetscReal)*nqp*nen*dim*dim);CHKERRQ(ierr);
    ierr = PetscMemzero(element->shape[3],sizeof(PetscReal)*nqp*nen*dim*dim*dim);CHKERRQ(ierr);

    ierr = PetscMemzero(element->normal,sizeof(PetscReal)*nqp*dim);CHKERRQ(ierr);
  }
  { /* */
    PetscInt nen = element->nen;
    PetscInt dof = element->dof;
    element->nfix = 0;
    ierr = PetscMalloc1(nen*dof,PetscInt,   &element->ifix);CHKERRQ(ierr);
    ierr = PetscMalloc1(nen*dof,PetscScalar,&element->vfix);CHKERRQ(ierr);
    ierr = PetscMalloc1(nen*dof,PetscScalar,&element->ufix);CHKERRQ(ierr);
    element->nflux = 0;
    ierr = PetscMalloc1(nen*dof,PetscInt,   &element->iflux);CHKERRQ(ierr);
    ierr = PetscMalloc1(nen*dof,PetscScalar,&element->vflux);CHKERRQ(ierr);
  }
  ierr = IGAPointInit(element->iterator,element);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAGetElement"
PetscErrorCode IGAGetElement(IGA iga,IGAElement *element)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidPointer(element,2);
  *element = iga->iterator;
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGABeginElement"
PetscErrorCode IGABeginElement(IGA iga,IGAElement *element)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidPointer(element,2);
  IGACheckSetUp(iga,1);
  *element = iga->iterator;
  (*element)->index = -1;
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGANextElement"
PetscBool IGANextElement(IGA iga,IGAElement element)
{
  PetscInt i,dim  = element->dim;
  PetscInt *start = element->start;
  PetscInt *width = element->width;
  PetscInt *ID    = element->ID;
  PetscInt index,coord;
  PetscErrorCode ierr;
  /* */
  element->nval = 0;
  element->nvec = 0;
  element->nmat = 0;
  element->boundary_id = -1;
  element->interface_id = -1;
  /* */
  index = ++element->index;
  if (PetscUnlikely(index >= element->count)) {
    element->index = -1;
    return PETSC_FALSE;
  }
  for (i=0; i<dim; i++) {
    coord = index % width[i];
    index = (index - coord) / width[i];
    ID[i] = coord + start[i];
  }
  /* */
#undef  CHKERRRETURN
#define CHKERRRETURN(n,r) do{if(PetscUnlikely(n)){CHKERRCONTINUE(n);return(r);}}while(0)
  ierr = IGAElementBuildMapping(element);  CHKERRRETURN(ierr,PETSC_FALSE);
  ierr = IGAElementBuildGeometry(element); CHKERRRETURN(ierr,PETSC_FALSE);
  ierr = IGAElementBuildProperty(element); CHKERRRETURN(ierr,PETSC_FALSE);
  ierr = IGAElementBuildFix(element);      CHKERRRETURN(ierr,PETSC_FALSE);
#undef  CHKERRRETURN
  return PETSC_TRUE;
}

#undef  __FUNCT__
#define __FUNCT__ "IGAEndElement"
PetscErrorCode IGAEndElement(IGA iga,IGAElement *element)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidPointer(element,2);
  PetscValidPointer(*element,2);
  if (PetscUnlikely((*element)->index != -1)) {
    (*element)->index = -1;
    *element = NULL;
    PetscFunctionReturn(PETSC_ERR_PLIB);
  }
  *element = NULL;
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAElementNextUserOps"
PetscBool IGAElementNextUserOps(IGAElement element,IGAUserOps *userops)
{
  IGA      iga = element->parent;
  PetscInt dim = element->dim;
  while (++element->boundary_id < 2*dim) {
    PetscInt i = element->boundary_id / 2;
    PetscInt s = element->boundary_id % 2;
    PetscInt e = s ? element->sizes[i]-1 : 0;
    if (element->ID[i] != e) continue;
    if (!iga->boundary[i][s]->userops) continue;
    *userops = iga->boundary[i][s]->userops;
    element->atboundary = PETSC_TRUE;
    return PETSC_TRUE;
  }
  while (++element->interface_id < 2*dim) {
//	PetscPrintf(PETSC_COMM_WORLD,"inside \n");
    PetscInt Nelt=3;
	PetscInt i = element->interface_id / 2;
    PetscInt s = element->interface_id % 2;
//    PetscInt e = s ? element->sizes[i]-1 : 0;
    if (iga->dim ==2) if (element->ID[i] != element->sizes[i]-Nelt || element->ID[0] > element->sizes[0]-Nelt || element->ID[1] > element->sizes[1]-Nelt) continue;
    if (iga->dim ==3) if (element->ID[i] != element->sizes[i]-Nelt || element->ID[0] > element->sizes[0]-Nelt || element->ID[1] > element->sizes[1]-Nelt || element->ID[2] > element->sizes[2]-Nelt) continue;
    if (!iga->interface[i][s]->userops) continue;
//    PetscPrintf(PETSC_COMM_WORLD,"mega inside \n");
//    PetscPrintf(PETSC_COMM_WORLD,"id:%d \n",element->ID[i]);
//    PetscPrintf(PETSC_COMM_WORLD,"size -Nelt:%d \n",element->sizes[i]-3);
//    	  	PetscPrintf(PETSC_COMM_WORLD,"ii:%d \n",i);
//    	  	PetscPrintf(PETSC_COMM_WORLD,"ss:%d \n",s);

    *userops = iga->interface[i][s]->userops;
    element->atinterface = PETSC_TRUE;
    return PETSC_TRUE;
  }
  if (element->boundary_id++ == 2*dim || element->interface_id++ == 2*dim) {
    *userops = iga->userops;
    element->atboundary = PETSC_FALSE;
    element->atinterface = PETSC_FALSE;
    return PETSC_TRUE;
  }

  *userops = 0;
  element->atboundary  = PETSC_FALSE;
  element->boundary_id = -1;
  element->atinterface = PETSC_FALSE;
  element->interface_id = -1;
  return PETSC_FALSE;
}

#undef  __FUNCT__
#define __FUNCT__ "IGAElementGetPoint"
PetscErrorCode IGAElementGetPoint(IGAElement element,IGAPoint *point)
{
  PetscFunctionBegin;
  PetscValidPointer(element,1);
  PetscValidPointer(point,2);
  *point = element->iterator;
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAElementBeginPoint"
PetscErrorCode IGAElementBeginPoint(IGAElement element,IGAPoint *point)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidPointer(element,1);
  PetscValidPointer(point,2);
  if (PetscUnlikely(element->index < 0))
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Must call during element loop");

  *point = element->iterator;
  (*point)->index = -1;
  (*point)->count = element->nqp;

  (*point)->neq = element->neq;
  (*point)->nen = element->nen;
  (*point)->dof = element->dof;
  (*point)->dim = element->dim;
  (*point)->nsd = element->nsd;
  (*point)->npd = element->npd;

  if (PetscLikely(!element->atboundary) && PetscLikely(!element->atinterface)) {
//    PetscPrintf(PETSC_COMM_WORLD,"El:%d \n",element->index);
    ierr = IGAElementBuildQuadrature(element);CHKERRQ(ierr);
    ierr = IGAElementBuildShapeFuns(element);CHKERRQ(ierr);
  } else if (element->atboundary) {
//	PetscPrintf(PETSC_COMM_WORLD,"El Boundary:%d \n",element->index);
    PetscInt i = element->boundary_id / 2;      //direction
    PetscInt s = element->boundary_id % 2;      //side
//	PetscPrintf(PETSC_COMM_WORLD,"i:%d \n",i);
//	PetscPrintf(PETSC_COMM_WORLD,"s:%d \n",s);
    (*point)->count = element->nqp / element->BD[i]->nqp;
    ierr = IGAElementBuildQuadratureAtBoundary(element,i,s);CHKERRQ(ierr);
    ierr = IGAElementBuildShapeFunsAtBoundary (element,i,s);CHKERRQ(ierr);
 } else if (element->atinterface){
//	  PetscPrintf(PETSC_COMM_WORLD,"El Interface:%d \n",element->index);
	PetscInt i = element->interface_id / 2;      //direction
	PetscInt s = element->interface_id % 2;      //side
//	  	PetscPrintf(PETSC_COMM_WORLD,"i:%d \n",i);
//	  	PetscPrintf(PETSC_COMM_WORLD,"s:%d \n",s);
	(*point)->count = element->nqp / element->BD[i]->nqp;
	ierr = IGAElementBuildQuadratureAtBoundary(element,i,s);CHKERRQ(ierr);
	ierr = IGAElementBuildShapeFunsAtBoundary (element,i,s);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAElementNextPoint"
PetscBool IGAElementNextPoint(IGAElement element,IGAPoint point)
{
  PetscInt nen = point->nen;
  PetscInt dim = point->dim;
  PetscInt nsd = point->nsd;
  PetscInt index;
  /* */
  point->nvec = 0;
  point->nmat = 0;
  /* */
  index = ++point->index;
  if (PetscUnlikely(index == 0))            goto start;
  if (PetscUnlikely(index >= point->count)) goto stop;

  point->weight   += 1;
  point->detJac   += 1;

  point->point    += dim;
  point->basis[0] += nen;
  point->basis[1] += nen*dim;
  point->basis[2] += nen*dim*dim;
  point->basis[3] += nen*dim*dim*dim;

  point->detX     += 1;
  point->gradX[0] += dim*dim;
  point->gradX[1] += dim*dim;
  point->shape[0] += nen;
  point->shape[1] += nen*dim;
  point->shape[2] += nen*dim*dim;
  point->shape[3] += nen*dim*dim*dim;

  point->normal   += dim;

  return PETSC_TRUE;

 start:

  point->geometry    = element->geometryX;
  point->geometryINI = element->geometryINI;
  point->rational    = element->rationalW;
  point->property    = element->propertyA;
  if (!element->geometry)
    point->geometryINI = NULL;
  if (!element->geometry)
    point->geometry = NULL;
  if (!element->rational)
    point->rational = NULL;
  if (!element->property)
    point->property = NULL;

  point->weight   = element->weight;
  point->detJac   = element->detJac;

  point->point    = element->point;
  point->basis[0] = element->basis[0];
  point->basis[1] = element->basis[1];
  point->basis[2] = element->basis[2];
  point->basis[3] = element->basis[3];

  if (element->geometry && dim == nsd) { /* XXX */
    point->detX     = element->detX;
    point->gradX[0] = element->gradX[0];
    point->gradX[1] = element->gradX[1];
    point->shape[0] = element->shape[0];
    point->shape[1] = element->shape[1];
    point->shape[2] = element->shape[2];
    point->shape[3] = element->shape[3];
  } else {
    point->shape[0] = element->basis[0];
    point->shape[1] = element->basis[1];
    point->shape[2] = element->basis[2];
    point->shape[3] = element->basis[3];
  }

  point->normal = element->normal;

  return PETSC_TRUE;

 stop:

  point->index = -1;
  return PETSC_FALSE;

}

#undef  __FUNCT__
#define __FUNCT__ "IGAElementEndPoint"
PetscErrorCode IGAElementEndPoint(IGAElement element,IGAPoint *point)
{
  PetscFunctionBegin;
  PetscValidPointer(element,1);
  PetscValidPointer(point,2);
  PetscValidPointer(*point,2);
  if (PetscUnlikely((*point)->index != -1)) {
    (*point)->index = -1;
    PetscFunctionReturn(PETSC_ERR_PLIB);
  }
  *point = NULL;
  /* XXX */
  if (PetscLikely(!element->collocation)) PetscFunctionReturn(0);
  if (PetscLikely(!element->atboundary))  PetscFunctionReturn(0);
  if (PetscLikely(!element->atinterface))  PetscFunctionReturn(0);
  element->atboundary  = PETSC_FALSE;
  element->atinterface = PETSC_FALSE;
  element->boundary_id = 2*element->dim;
  element->interface_id = 2*element->dim;
  /* XXX */
 PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAElementGetParent"
PetscErrorCode IGAElementGetParent(IGAElement element,IGA *parent)
{
  PetscFunctionBegin;
  PetscValidPointer(element,1);
  PetscValidIntPointer(parent,2);
  *parent = element->parent;
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAElementGetIndex"
PetscErrorCode IGAElementGetIndex(IGAElement element,PetscInt *index)
{
  PetscFunctionBegin;
  PetscValidPointer(element,1);
  PetscValidIntPointer(index,2);
  *index = element->index;
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAElementGetCount"
PetscErrorCode IGAElementGetCount(IGAElement element,PetscInt *count)
{
  PetscFunctionBegin;
  PetscValidPointer(element,1);
  PetscValidIntPointer(count,2);
  *count = element->count;
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAElementGetSizes"
PetscErrorCode IGAElementGetSizes(IGAElement element,PetscInt *neq,PetscInt *nen,PetscInt *dof)
{
  PetscFunctionBegin;
  PetscValidPointer(element,1);
  if (neq) PetscValidIntPointer(neq,2);
  if (nen) PetscValidIntPointer(nen,3);
  if (dof) PetscValidIntPointer(dof,4);
  if (neq) *neq = element->neq;
  if (nen) *nen = element->nen;
  if (dof) *dof = element->dof;
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAElementGetMapping"
PetscErrorCode IGAElementGetMapping(IGAElement element,PetscInt *nen,const PetscInt *mapping[])
{
  PetscFunctionBegin;
  PetscValidPointer(element,1);
  if (nen)     PetscValidIntPointer(nen,3);
  if (mapping) PetscValidPointer(mapping,3);
  if (nen)     *nen     = element->nen;
  if (mapping) *mapping = element->mapping;
  PetscFunctionReturn(0);
}

/* -------------------------------------------------------------------------- */

#undef  __FUNCT__
#define __FUNCT__ "IGAElementBuildMapping"
PetscErrorCode IGAElementBuildMapping(IGAElement element)
{
  PetscFunctionBegin;
  PetscValidPointer(element,1);
  if (PetscUnlikely(element->index < 0))
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Must call during element loop");
  { /* */
    IGA      iga = element->parent;
    IGABasis *BD = element->BD;
    PetscInt *ID = element->ID;
    PetscInt ia, inen = BD[0]->nen, ioffset = BD[0]->offset[ID[0]];
    PetscInt ja, jnen = BD[1]->nen, joffset = BD[1]->offset[ID[1]];
    PetscInt ka, knen = BD[2]->nen, koffset = BD[2]->offset[ID[2]];
    PetscInt *start = iga->node_gstart, *width = iga->node_gwidth;
    PetscInt istart = start[0]/*istride = 1*/;
    PetscInt jstart = start[1], jstride = width[0];
    PetscInt kstart = start[2], kstride = width[0]*width[1];
    PetscInt a=0, *mapping = element->mapping;
    for (ka=0; ka<knen; ka++) {
      for (ja=0; ja<jnen; ja++) {
        for (ia=0; ia<inen; ia++) {
          PetscInt iA = (ioffset + ia) - istart;
          PetscInt jA = (joffset + ja) - jstart;
          PetscInt kA = (koffset + ka) - kstart;
          mapping[a++] = iA + jA*jstride + kA*kstride;  //Relates each node of the element to the "local" processor node
        }
      }
    }
    if (PetscUnlikely(element->collocation)) {
      PetscInt iA = ID[0] - istart;
      PetscInt jA = ID[1] - jstart;
      PetscInt kA = ID[2] - kstart;
      element->rowmap[0] = iA + jA*jstride + kA*kstride;
    }
  }
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAElementBuildGeometry"
PetscErrorCode IGAElementBuildGeometry(IGAElement element)
{
  IGA iga;
  PetscFunctionBegin;
  PetscValidPointer(element,1);
  if (PetscUnlikely(element->index < 0))
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Must call during element loop");
  iga = element->parent;
  element->geometry = iga->geometry ? PETSC_TRUE : PETSC_FALSE;
  element->rational = iga->rational ? PETSC_TRUE : PETSC_FALSE;
  if (element->geometry || element->rational) {
    const PetscInt  *map = element->mapping;
    const PetscReal *arrayX   = iga->geometryX;
    const PetscReal *arrayINI = iga->geometryINI;
    const PetscReal *arrayW   = iga->rationalW;
    PetscReal *X   = element->geometryX;
    PetscReal *INI = element->geometryINI;
    PetscReal *W   = element->rationalW;
    PetscInt a,nen = element->nen;
    PetscInt i,nsd = element->nsd;
    if (element->geometry)
      for (a=0; a<nen; a++)
        for (i=0; i<nsd; i++) {
          X[i + a*nsd] = arrayX[i + map[a]*nsd];
          INI[i + a*nsd] = arrayINI[i + map[a]*nsd];
        }
    if (element->rational)
      for (a=0; a<nen; a++)
        W[a] = arrayW[map[a]];
  }
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAElementBuildProperty"
PetscErrorCode IGAElementBuildProperty(IGAElement element)
{
  IGA iga;
  PetscFunctionBegin;
  PetscValidPointer(element,1);
  if (PetscUnlikely(element->index < 0))
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Must call during element loop");
  iga = element->parent;
  element->property = iga->property ? PETSC_TRUE: PETSC_FALSE;
  if (element->property) {
    const PetscInt *map = element->mapping;
    const PetscScalar *arrayA = iga->propertyA;
    PetscScalar *A = element->propertyA;
    PetscInt a,nen = element->nen;
    PetscInt i,npd = element->npd;
    for (a=0; a<nen; a++)
      for (i=0; i<npd; i++)
        A[i + a*npd] = arrayA[i + map[a]*npd];
  }
  PetscFunctionReturn(0);
}

/* -------------------------------------------------------------------------- */

EXTERN_C_BEGIN
extern void IGA_Quadrature_1D(PetscInt,const PetscReal[],const PetscReal[],const PetscReal*,
                              PetscReal[],PetscReal[],PetscReal[]);
extern void IGA_Quadrature_2D(PetscInt,const PetscReal[],const PetscReal[],const PetscReal*,
                              PetscInt,const PetscReal[],const PetscReal[],const PetscReal*,
                              PetscReal[],PetscReal[],PetscReal[]);
extern void IGA_Quadrature_3D(PetscInt,const PetscReal[],const PetscReal[],const PetscReal*,
                              PetscInt,const PetscReal[],const PetscReal[],const PetscReal*,
                              PetscInt,const PetscReal[],const PetscReal[],const PetscReal*,
                              PetscReal[],PetscReal[],PetscReal[]);
EXTERN_C_END

EXTERN_C_BEGIN
extern void IGA_BasisFuns_1D(PetscInt,PetscInt,const PetscReal[],
                             PetscInt,PetscInt,PetscInt,const PetscReal[],
                             PetscReal[],PetscReal[],PetscReal[],PetscReal[]);
extern void IGA_BasisFuns_2D(PetscInt,PetscInt,const PetscReal[],
                             PetscInt,PetscInt,PetscInt,const PetscReal[],
                             PetscInt,PetscInt,PetscInt,const PetscReal[],
                             PetscReal[],PetscReal[],PetscReal[],PetscReal[]);
extern void IGA_BasisFuns_3D(PetscInt,PetscInt,const PetscReal[],
                             PetscInt,PetscInt,PetscInt,const PetscReal[],
                             PetscInt,PetscInt,PetscInt,const PetscReal[],
                             PetscInt,PetscInt,PetscInt,const PetscReal[],
                             PetscReal[],PetscReal[],PetscReal[],PetscReal[]);
EXTERN_C_END

EXTERN_C_BEGIN
extern void IGA_ShapeFuns_1D(PetscInt,PetscInt,PetscInt,const PetscReal[],
                             const PetscReal[],const PetscReal[],const PetscReal[],const PetscReal[],
                             PetscReal[],PetscReal[],PetscReal[],PetscReal[],
                             PetscReal[],PetscReal[],PetscReal[]);
extern void IGA_ShapeFuns_2D(PetscInt,PetscInt,PetscInt,const PetscReal[],
                             const PetscReal[],const PetscReal[],const PetscReal[],const PetscReal[],
                             PetscReal[],PetscReal[],PetscReal[],PetscReal[],
                             PetscReal[],PetscReal[],PetscReal[]);
extern void IGA_ShapeFuns_3D(PetscInt,PetscInt,PetscInt,const PetscReal[],
                             const PetscReal[],const PetscReal[],const PetscReal[],const PetscReal[],
                             PetscReal[],PetscReal[],PetscReal[],PetscReal[],
                             PetscReal[],PetscReal[],PetscReal[]);
EXTERN_C_END

#define IGA_Quadrature_ARGS(ID,BD,i) \
  BD[i]->nqp,BD[i]->point+ID[i]*BD[i]->nqp,BD[i]->weight,BD[i]->detJ+ID[i]

#define IGA_BasisFuns_ARGS(ID,BD,i) \
  BD[i]->nqp,BD[i]->nen,BD[i]->d,BD[i]->value+ID[i]*BD[i]->nqp*BD[i]->nen*(BD[i]->d+1)

#undef  __FUNCT__
#define __FUNCT__ "IGAElementBuildQuadrature"
PetscErrorCode IGAElementBuildQuadrature(IGAElement element)
{
  PetscFunctionBegin;
  PetscValidPointer(element,1);
  if (PetscUnlikely(element->index < 0))
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Must call during element loop");
  {
    IGABasis *BD = element->BD;
    PetscInt *ID = element->ID;
    PetscReal *u = element->point;
    PetscReal *w = element->weight;
    PetscReal *J = element->detJac;
    switch (element->dim) {
    case 3: IGA_Quadrature_3D(IGA_Quadrature_ARGS(ID,BD,0),
                              IGA_Quadrature_ARGS(ID,BD,1),
                              IGA_Quadrature_ARGS(ID,BD,2),
                              u,w,J); break;
    case 2: IGA_Quadrature_2D(IGA_Quadrature_ARGS(ID,BD,0),
                              IGA_Quadrature_ARGS(ID,BD,1),
                              u,w,J); break;
    case 1: IGA_Quadrature_1D(IGA_Quadrature_ARGS(ID,BD,0),
                              u,w,J); break;
    }
  }
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAElementBuildShapeFuns"
PetscErrorCode IGAElementBuildShapeFuns(IGAElement element)
{
  PetscFunctionBegin;
  PetscValidPointer(element,1);
  if (PetscUnlikely(element->index < 0))
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Must call during element loop");
  {
    IGABasis  *BD = element->BD;
    PetscInt  *ID = element->ID;
    PetscInt  ord = element->parent->order;
    PetscInt  rat = element->rational;
    PetscReal *W  = element->rationalW;
    PetscReal **N = element->basis;
    switch (element->dim) {
    case 3: IGA_BasisFuns_3D(ord,rat,W,
                             IGA_BasisFuns_ARGS(ID,BD,0),
                             IGA_BasisFuns_ARGS(ID,BD,1),
                             IGA_BasisFuns_ARGS(ID,BD,2),
                             N[0],N[1],N[2],N[3]); break;
    case 2: IGA_BasisFuns_2D(ord,rat,W,
                             IGA_BasisFuns_ARGS(ID,BD,0),
                             IGA_BasisFuns_ARGS(ID,BD,1),
                             N[0],N[1],N[2],N[3]); break;
    case 1: IGA_BasisFuns_1D(ord,rat,W,
                             IGA_BasisFuns_ARGS(ID,BD,0),
                             N[0],N[1],N[2],N[3]); break;
    }
    {
      PetscInt nqp = element->nqp;
      PetscInt dim = element->dim;
      PetscReal *n = element->normal;
      (void)PetscMemzero(n,nqp*dim*sizeof(PetscReal));
    }
  }

  if (element->dim == element->nsd) /* XXX */
  if (element->geometry) {
    PetscInt q;
    PetscInt ord  = element->parent->order;
    PetscInt nqp  = element->nqp;

    PetscInt nen  = element->nen;


//    PetscReal *X;
//    if (element->parent->cruz == 1){
//    	X  = element->geometryX;
//    } else{
//    	if (element->parent->type[element->index] == 1 ) {  //Solid
//    	   	X  = element->geometryINI;
//    	} else {											//Fluid
//        	X  = element->geometryX;
//    	}
//    }

    PetscReal *X  = element->geometryX;
//    PetscReal *INI= element->geometryINI;
//    PetscReal *X  = element->geometryINI;

    PetscReal **M = element->basis;
    PetscReal **N = element->shape;
    PetscReal *J  = element->detX;
    PetscReal *F  = element->gradX[0];
    PetscReal *G  = element->gradX[1];
    switch (element->dim) {
    case 3: IGA_ShapeFuns_3D(ord,nqp,nen,X,
                             M[0],M[1],M[2],M[3],
                             N[0],N[1],N[2],N[3],
                             J,F,G); break;
    case 2: IGA_ShapeFuns_2D(ord,nqp,nen,X,
                             M[0],M[1],M[2],M[3],
                             N[0],N[1],N[2],N[3],
                             J,F,G); break;
    case 1: IGA_ShapeFuns_1D(ord,nqp,nen,X,
                             M[0],M[1],M[2],M[3],
                             N[0],N[1],N[2],N[3],
                             J,F,G); break;
    }
    for (q=0; q<nqp; q++)
      element->detJac[q] *= J[q];
  }
  PetscFunctionReturn(0);
}

#define IGA_Quadrature_BNDR(ID,BD,i,s) \
  1,&BD[i]->bnd_point[s],&BD[i]->bnd_weight[s],&BD[i]->bnd_detJ[s]

#define IGA_BasisFuns_BNDR(ID,BD,i,s) \
  1,BD[i]->nen,BD[i]->d,BD[i]->bnd_value[s]

EXTERN_C_BEGIN
extern void IGA_GetNormal(PetscInt dim,PetscInt dir,PetscInt side,const PetscReal F[],PetscReal *dS,PetscReal N[]);
EXTERN_C_END

#undef  __FUNCT__
#define __FUNCT__ "IGAElementBuildQuadratureAtBoundary"
PetscErrorCode IGAElementBuildQuadratureAtBoundary(IGAElement element,PetscInt dir,PetscInt side)
{
  PetscFunctionBegin;
  PetscValidPointer(element,1);
  if (PetscUnlikely(element->index < 0))
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Must call during element loop");
  {
    IGABasis *BD = element->BD;
    PetscInt *ID = element->ID;
    PetscReal *u = element->point;
    PetscReal *w = element->weight;
    PetscReal *J = element->detJac;
    switch (element->dim) {
    case 3:
      switch (dir) {
      case 0: IGA_Quadrature_3D(IGA_Quadrature_BNDR(ID,BD,0,side),
                                IGA_Quadrature_ARGS(ID,BD,1),
                                IGA_Quadrature_ARGS(ID,BD,2),
                                u,w,J); break;
      case 1: IGA_Quadrature_3D(IGA_Quadrature_ARGS(ID,BD,0),
                                IGA_Quadrature_BNDR(ID,BD,1,side),
                                IGA_Quadrature_ARGS(ID,BD,2),
                                u,w,J); break;
      case 2: IGA_Quadrature_3D(IGA_Quadrature_ARGS(ID,BD,0),
                                IGA_Quadrature_ARGS(ID,BD,1),
                                IGA_Quadrature_BNDR(ID,BD,2,side),
                                u,w,J); break;
      } break;
    case 2:
      switch (dir) {
      case 0: IGA_Quadrature_2D(IGA_Quadrature_BNDR(ID,BD,0,side),
                                IGA_Quadrature_ARGS(ID,BD,1),
                                u,w,J); break;
      case 1: IGA_Quadrature_2D(IGA_Quadrature_ARGS(ID,BD,0),
                                IGA_Quadrature_BNDR(ID,BD,1,side),
                                u,w,J); break;
      } break;
    case 1:
      switch (dir) {
      case 0: IGA_Quadrature_1D(IGA_Quadrature_BNDR(ID,BD,0,side),
                                u,w,J); break;
      } break;
    }
  }
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAElementBuildShapeFunsAtBoundary"
PetscErrorCode IGAElementBuildShapeFunsAtBoundary(IGAElement element,PetscInt dir,PetscInt side)
{
  PetscFunctionBegin;
  PetscValidPointer(element,1);
  if (PetscUnlikely(element->index < 0))
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Must call during element loop");
  {
    IGABasis  *BD = element->BD;
    PetscInt  *ID = element->ID;
    PetscInt  ord = element->parent->order;
    PetscInt  rat = element->rational;
    PetscReal *W  = element->rationalW;
    PetscReal **N = element->basis;
    switch (element->dim) {
    case 3:
      switch (dir) {
      case 0: IGA_BasisFuns_3D(ord,rat,W,
                               IGA_BasisFuns_BNDR(ID,BD,0,side),
                               IGA_BasisFuns_ARGS(ID,BD,1),
                               IGA_BasisFuns_ARGS(ID,BD,2),
                               N[0],N[1],N[2],N[3]); break;
      case 1: IGA_BasisFuns_3D(ord,rat,W,
                               IGA_BasisFuns_ARGS(ID,BD,0),
                               IGA_BasisFuns_BNDR(ID,BD,1,side),
                               IGA_BasisFuns_ARGS(ID,BD,2),
                               N[0],N[1],N[2],N[3]); break;
      case 2: IGA_BasisFuns_3D(ord,rat,W,
                               IGA_BasisFuns_ARGS(ID,BD,0),
                               IGA_BasisFuns_ARGS(ID,BD,1),
                               IGA_BasisFuns_BNDR(ID,BD,2,side),
                               N[0],N[1],N[2],N[3]); break;
      } break;
    case 2:
      switch (dir) {
      case 0: IGA_BasisFuns_2D(ord,rat,W,
                               IGA_BasisFuns_BNDR(ID,BD,0,side),
                               IGA_BasisFuns_ARGS(ID,BD,1),
                               N[0],N[1],N[2],N[3]); break;
      case 1: IGA_BasisFuns_2D(ord,rat,W,
                               IGA_BasisFuns_ARGS(ID,BD,0),
                               IGA_BasisFuns_BNDR(ID,BD,1,side),
                               N[0],N[1],N[2],N[3]); break;
      } break;
    case 1:
      switch (dir) {
      case 0: IGA_BasisFuns_1D(ord,rat,W,
                               IGA_BasisFuns_BNDR(ID,BD,0,side),
                               N[0],N[1],N[2],N[3]); break;
      } break;
    }
    {
      PetscInt q;
      PetscInt nqp = element->nqp;
      PetscInt dim = element->dim;
      PetscReal *n = element->normal;
      (void)PetscMemzero(n,nqp*dim*sizeof(PetscReal));
      for (q=0; q<nqp; q++)
        n[q*dim+dir] = side ? 1.0 : -1.0;
    }
  }

  if (element->dim == element->nsd) /* XXX */
  if (element->geometry) {
    PetscInt q;
    PetscInt dim  = element->dim;
    PetscInt ord  = element->parent->order;
    PetscInt nqp  = element->nqp;
    PetscInt nen  = element->nen;
/*
    PetscReal *X;
    if (element->parent->property) {
    	X  = element->geometryX;
    } else{
    	X  = element->geometryINI;

    }
*/

	PetscReal *X  = element->geometryX;
//    PetscReal *INI= element->geometryINI;
//    PetscReal *X= element->geometryINI;

    PetscReal **M = element->basis;
    PetscReal **N = element->shape;
    PetscReal *J  = element->detX;
    PetscReal *F  = element->gradX[0];
    PetscReal *G  = element->gradX[1];
    PetscReal *n  = element->normal;
    switch (dim) {
    case 3: IGA_ShapeFuns_3D(ord,nqp,nen,X,
                             M[0],M[1],M[2],M[3],
                             N[0],N[1],N[2],N[3],
                             J,F,G); break;
    case 2: IGA_ShapeFuns_2D(ord,nqp,nen,X,
                             M[0],M[1],M[2],M[3],
                             N[0],N[1],N[2],N[3],
                             J,F,G); break;
    case 1: IGA_ShapeFuns_1D(ord,nqp,nen,X,
                             M[0],M[1],M[2],M[3],
                             N[0],N[1],N[2],N[3],
                             J,F,G); break;
    }
    for (q=0; q<nqp; q++) {
      PetscReal dS;
      IGA_GetNormal(dim,dir,side,&F[q*dim*dim],&dS,&n[q*dim]);
      element->detJac[q] *= dS;
    }
/*	PetscPrintf(PETSC_COMM_WORLD,"Element: %d",element->index);
	PetscPrintf(PETSC_COMM_WORLD," normal: %e ",n[0]);
	PetscPrintf(PETSC_COMM_WORLD," normal: %e \n",n[1]);   */
  }
  PetscFunctionReturn(0);
}

/* -------------------------------------------------------------------------- */

#undef  __FUNCT__
#define __FUNCT__ "IGAElementGetWorkVal"
PetscErrorCode IGAElementGetWorkVal(IGAElement element,PetscScalar *U[])
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidPointer(element,1);
  PetscValidPointer(U,2);
  if (PetscUnlikely(element->index < 0))
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Must call during element loop");
  {
    size_t MAX_WORK_VAL = sizeof(element->wval)/sizeof(PetscScalar*);
    PetscInt n = element->nen * element->dof;
    if (PetscUnlikely(element->nval >= (PetscInt)MAX_WORK_VAL))
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Too many work values requested");
    if (PetscUnlikely(!element->wval[element->nval])) {
      ierr = PetscMalloc1(n,PetscScalar,&element->wval[element->nval]);CHKERRQ(ierr);
    }
    *U = element->wval[element->nval++];
    ierr = PetscMemzero(*U,n*sizeof(PetscScalar));CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAElementGetWorkVec"
PetscErrorCode IGAElementGetWorkVec(IGAElement element,PetscScalar *V[])
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidPointer(element,1);
  PetscValidPointer(V,2);
  if (PetscUnlikely(element->index < 0))
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Must call during element loop");
  {
    size_t MAX_WORK_VEC = sizeof(element->wvec)/sizeof(PetscScalar*);
    PetscInt m = element->neq * element->dof;
    PetscInt n = element->nen * element->dof;
    if (PetscUnlikely(element->nvec >= (PetscInt)MAX_WORK_VEC))
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Too many work vectors requested");
    if (PetscUnlikely(!element->wvec[element->nvec])) {
      ierr = PetscMalloc1(n,PetscScalar,&element->wvec[element->nvec]);CHKERRQ(ierr);
    }
    *V = element->wvec[element->nvec++];
    ierr = PetscMemzero(*V,m*sizeof(PetscScalar));CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAElementGetWorkMat"
PetscErrorCode IGAElementGetWorkMat(IGAElement element,PetscScalar *M[])
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidPointer(element,1);
  PetscValidPointer(M,2);
  if (PetscUnlikely(element->index < 0))
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Must call during element loop");
  {
    size_t MAX_WORK_MAT = sizeof(element->wmat)/sizeof(PetscScalar*);
    PetscInt m = element->neq * element->dof;
    PetscInt n = element->nen * element->dof;
    if (PetscUnlikely(element->nmat >= (PetscInt)MAX_WORK_MAT))
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Too many work matrices requested");
    if (PetscUnlikely(!element->wmat[element->nmat])) {
      ierr = PetscMalloc1(n*n,PetscScalar,&element->wmat[element->nmat]);CHKERRQ(ierr);
    }
    *M = element->wmat[element->nmat++];
    ierr = PetscMemzero(*M,m*n*sizeof(PetscScalar));CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/* -------------------------------------------------------------------------- */

#undef  __FUNCT__
#define __FUNCT__ "IGAElementGetValues"
PetscErrorCode IGAElementGetValues(IGAElement element,const PetscScalar arrayU[],PetscScalar U[])
{
  PetscFunctionBegin;
  PetscValidPointer(element,1);
  PetscValidScalarPointer(arrayU,2);
  PetscValidScalarPointer(U,3);
  {
    PetscInt nen = element->nen;
    PetscInt dof = element->dof;
    PetscInt *mapping = element->mapping;
    PetscInt a,i,pos=0;
    for (a=0; a<nen; a++) {
      const PetscScalar *u = arrayU + mapping[a]*dof;
      for (i=0; i<dof; i++)
        U[pos++] = u[i]; /* XXX Use PetscMemcpy() ?? */
    }
  }
  PetscFunctionReturn(0);
}

/* -------------------------------------------------------------------------- */

EXTERN_C_BEGIN
extern void IGA_BoundaryArea_2D(const PetscInt[],PetscInt,PetscInt,
                                PetscInt,const PetscReal[],
                                PetscInt,const PetscReal[],
                                PetscInt,const PetscReal[],PetscInt,PetscInt,const PetscReal[],
                                PetscReal*);
extern void IGA_BoundaryArea_3D(const PetscInt[],PetscInt,PetscInt,
                                PetscInt,const PetscReal[],
                                PetscInt,const PetscReal[],
                                PetscInt,const PetscReal[],PetscInt,PetscInt,const PetscReal[],
                                PetscInt,const PetscReal[],PetscInt,PetscInt,const PetscReal[],
                                PetscReal*);
EXTERN_C_END

static PetscReal BoundaryArea(IGAElement element,PetscInt dir,PetscInt side)
{
  PetscReal A = 1.0;
  PetscInt *ID = element->ID;
  IGABasis *BD = element->BD;
  PetscInt i,dim = element->dim;


  PetscPrintf(PETSC_COMM_WORLD,"ElementBoundaryArea:%d \n",element->index);


  if (dim == 1) return A;
  for (i=0; i<dim; i++)
    if (i != dir) {
      PetscReal L = BD[i]->detJ[ID[i]];
      PetscInt  n = BD[i]->nen;
      A *= L/n;
    }
  if (!element->geometry) {
    A *= (dim==2) ? 2 : 4; /* sum(W) = 2 */
  } else {
    PetscInt shape[3] = {1,1,1};
    PetscInt k,nqp[3],nen[3],ndr[3];
    PetscReal *W[3],*N[3],dS = 1.0;
    for (i=0; i<dim; i++)
      shape[i] = BD[i]->nen;
    for (k=0,i=0; i<dim; i++) {
      if (i == dir) continue;
      nqp[k] = BD[i]->nqp;
      nen[k] = BD[i]->nen;
      ndr[k] = BD[i]->d;
      W[k]   = BD[i]->weight;
      N[k]   = BD[i]->value+ID[i]*nqp[k]*nen[k]*(ndr[k]+1);
      k++;
    }
    switch (dim) {
    case 2: IGA_BoundaryArea_2D(shape,dir,side,
                                element->geometry,element->geometryX,
                                element->rational,element->rationalW,
                                nqp[0],W[0],nen[0],ndr[0],N[0],
                                &dS); break;
    case 3: IGA_BoundaryArea_3D(shape,dir,side,
                                element->geometry,element->geometryX,
                                element->rational,element->rationalW,
                                nqp[0],W[0],nen[0],ndr[0],N[0],
                                nqp[1],W[1],nen[1],ndr[1],N[1],
                                &dS);break;
    }
    A *= dS;
  }
  return A;
}

static void AddFixa(IGAElement element,IGABoundary b,PetscInt a)
{
  if (b->count) {
    PetscInt dof = element->dof;
    PetscInt count = element->nfix;
    PetscInt *index = element->ifix;
    PetscScalar *value = element->vfix;
    PetscInt j,k,n = b->count;
    for (k=0; k<n; k++) {
      PetscInt idx = a*dof + b->field[k];
      PetscScalar val = b->value[k];
      for (j=0; j<count; j++)
        if (index[j] == idx) break;
      if (j == count) count++;
      index[j] = idx;
      value[j] = val;
    }
    element->nfix = count;
  }
}

static void AddFlux(IGAElement element,IGABoundary b,PetscInt a,PetscReal A)
{
  if (b->nload) {
    PetscInt dof = element->dof;
    PetscInt count = element->nflux;
    PetscInt *index = element->iflux;
    PetscScalar *value = element->vflux;
    PetscInt j,k,n = b->nload;
    for (k=0; k<n; k++) {
      PetscInt idx = a*dof + b->iload[k];
      PetscScalar val = b->vload[k];
      for (j=0; j<count; j++)
        if (index[j] == idx) break;
      if (j == count) value[count++] = 0.0;
      index[j] = idx;
      value[j]+= val*A;
    }
    element->nflux = count;
  }
}

static void BuildFix(IGAElement element,PetscInt dir,PetscInt side)
{
  IGABoundary b = element->parent->boundary[dir][side];
  if (b->count || b->nload) {
    PetscReal Area = b->nload ? BoundaryArea(element,dir,side) : 1.0;
    IGABasis *BD = element->BD;
    PetscInt S[3]={0,0,0},E[3]={1,1,1};
    PetscInt ia,ja,ka,jstride,kstride,a;
    PetscInt i,dim = element->dim;
    for (i=0; i<dim; i++) E[i] = BD[i]->nen;
    jstride = E[0]; kstride = E[0]*E[1];
    if (side) S[dir] = E[dir]-1;
    else      E[dir] = S[dir]+1;
    for (ka=S[2]; ka<E[2]; ka++)
      for (ja=S[1]; ja<E[1]; ja++)
        for (ia=S[0]; ia<E[0]; ia++)
          {
            a = ia + ja*jstride + ka*kstride;
            AddFixa(element,b,a);
            AddFlux(element,b,a,Area);
          }
  }
}

PETSC_STATIC_INLINE
IGABoundary AtBoundary(IGAElement element,PetscInt dir,PetscInt side)
{
  IGABoundary b = element->parent->boundary[dir][side];
  PetscInt e = side ? element->sizes[dir]-1 : 0;
  return (element->ID[dir] == e) ? b : NULL;
}

PETSC_STATIC_INLINE
PetscReal DOT(PetscInt dim,const PetscReal a[],const PetscReal b[])
{
  PetscInt i; PetscReal s = 0.0;
  for (i=0; i<dim; i++) s += a[i]*b[i];
  return s;
}

#undef  __FUNCT__
#define __FUNCT__ "IGAElementBuildFix"
PetscErrorCode IGAElementBuildFix(IGAElement element)
{
  PetscFunctionBegin;
  PetscValidPointer(element,1);
  if (PetscUnlikely(element->index < 0))
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Must call during element loop");
  if (PetscUnlikely(element->collocation)) goto collocation;
  element->nfix  = 0;
  element->nflux = 0;
  {
    IGAAxis  *AX = element->parent->axis;
    PetscInt *ID = element->ID;
    PetscInt i,dim = element->dim;
    for (i=0; i<dim; i++) {
      PetscBool w = AX[i]->periodic;
      PetscInt  e = element->sizes[i]-1; /* last element */
      if (ID[i] == 0 && !w) BuildFix(element,i,0);
      if (ID[i] == e && !w) BuildFix(element,i,1);
    }
  }
  PetscFunctionReturn(0);
 collocation:
  element->nfix  = 0;
  element->nflux = 0;
  {
    PetscInt L[3] = {PETSC_MIN_INT,PETSC_MIN_INT,PETSC_MIN_INT};
    PetscInt R[3] = {PETSC_MAX_INT,PETSC_MAX_INT,PETSC_MAX_INT};
    {
      IGAAxis  *AX = element->parent->axis;
      PetscInt i,dim = element->dim;
      for (i=0; i<dim; i++) {
        PetscBool w = AX[i]->periodic;
        PetscInt  n = element->sizes[i]-1; /* last node */
        L[i] = 0; if (!w) R[i] = n;
      }
    }
    {
      IGABoundary (*b)[2] = element->parent->boundary;
      IGABasis *BD = element->BD;
      PetscInt *ID = element->ID;
      PetscInt ia, inen = BD[0]->nen, ioffset = BD[0]->offset[ID[0]];
      PetscInt ja, jnen = BD[1]->nen, joffset = BD[1]->offset[ID[1]];
      PetscInt ka, knen = BD[2]->nen, koffset = BD[2]->offset[ID[2]];
      PetscInt a = 0;
      for (ka=0; ka<knen; ka++)
        for (ja=0; ja<jnen; ja++)
          for (ia=0; ia<inen; ia++)
            {
              PetscInt iA = ioffset + ia;
              PetscInt jA = joffset + ja;
              PetscInt kA = koffset + ka;
              /**/ if (iA == L[0]) AddFixa(element,b[0][0],a);
              else if (iA == R[0]) AddFixa(element,b[0][1],a);
              /**/ if (jA == L[1]) AddFixa(element,b[1][0],a);
              else if (jA == R[1]) AddFixa(element,b[1][1],a);
              /**/ if (kA == L[2]) AddFixa(element,b[2][0],a);
              else if (kA == R[2]) AddFixa(element,b[2][1],a);
              a++;
            }
    }
  }
  PetscFunctionReturn(0);
}



#undef  __FUNCT__
#define __FUNCT__ "IGAElementFixValues0_JB"
PetscErrorCode IGAElementFixValues0_JB(IGA iga, IGAElement element)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidPointer(element,2);
  {

	    PetscInt a,b;
	    PetscInt nen = element->nen;
	    PetscInt dim = element->dim;
	    PetscInt dof = element->dof;
	    PetscInt npd = element->npd;
//	    const PetscInt *map   = element->mapping;
//	    PetscInt       count  = element->nfix;
	    PetscInt       *index = element->ifix;
	    PetscScalar    *value = element->vfix;
	    PetscInt j=0.;


//	    PetscPrintf(PETSC_COMM_WORLD,"Element:%d \n",element->index);

/*	    PetscInt i;
	    for (i=0; i<nen; i++) {
	    for (b=0; b<npd; b++) {
	    	PetscPrintf(PETSC_COMM_WORLD,"Node, %d ELType: %e \n",i,element->propertyA[i*npd+b]);
        }
	    }
*/

	    if(iga->cruz == 1){



	//Displacements on the solid to the fluid mesh
	    for (a=0; a<nen; a++) {
	    for (b=0; b<dim; b++) {
	      index[j] = a*dim + b;
			if(element->propertyA[a*npd + 1] == 8.) {
				value[j] = element->propertyA[a*npd + 2 + b];
//		        PetscPrintf(PETSC_COMM_WORLD,"%d Solid: %e Indx: %d \n",j, value[j], index[j]);
				j++;
			}
	    }
	    }

	    // Symmetry conditions (comment for drop case)
	    for (b=0; b<dim; b++) {
	    for (a=0; a<nen; a++) {
			if(element->propertyA[a*npd + 1] != 8.) {
				if( element->geometryX[a*dim + b ] <= 1e-10) {
					index[j] = a*dim + b;
					value[j] = 0;
					j++;
			    }
			}

	    }
	    }


	    } else {

        //BC for temperature on the solid nodes (comment for drop case)
	    for (a=0; a<nen; a++) {
	      index[j] = a*dof + dim+2;
			if(element->propertyA[a*npd + 1] == 8.) {
//				value[j] = 0.;
				value[j] = iga->T_wall;  //323.5;
//		        PetscPrintf(PETSC_COMM_WORLD,"%d Val: %e Indx: %d \n",j, value[j], index[j]);
				j++;
			}
	    }


        //Symmetry conditions (comment for drop case)
	    for (b=0; b<dim; b++) {
	    for (a=0; a<nen; a++) {
				if( element->geometryX[a*dim + b ] <= 1e-10) {
					index[j] = a*dof + 1 + b;
					value[j] = 0;
					j++;
			    }

				//////////////////////////////// nsk_inter (v = 0 at the solid boundary)
//				if( element->geometryX[a*dim + b ] == 1.) {
//									index[j] = a*dof + 1 + b;
//									value[j] = 0;
//									j++;
//				}
				//////////////////////////////////////

	    }
	    }


	    }





	    element->nfix = j;
//	    PetscPrintf(PETSC_COMM_WORLD," Nfix:%d \n",element->nfix);

    }

  PetscFunctionReturn(0);

}


#undef  __FUNCT__
#define __FUNCT__ "IGAElementFixValues_JB"
PetscErrorCode IGAElementFixValues_JB(IGAElement element,PetscScalar U[],PetscScalar V[],PetscScalar A[])
{
  PetscFunctionBegin;
  PetscValidPointer(element,1);
  PetscValidScalarPointer(U,2);
  PetscValidScalarPointer(V,3);
  PetscValidScalarPointer(A,4);
  {
    PetscInt f,n;
    n = element->nfix;

    if(element->parent->cruz == 1) {


    for (f=0; f<n; f++) {
      PetscInt k = element->ifix[f];
      element->ufix[f] = U[k];
      U[k] = element->vfix[f]*0.666666;
//      V[k] = element->vfix[f]*0.666666*0.6666666/(0.340277*1e-2);       //Not used
//      A[k] = element->vfix[f]*0.833333/(0.340277*1e-4);                 //Not used

/*
 //#############################################
   PetscPrintf(PETSC_COMM_WORLD,"ifix:%d \n",element->ifix[f]);
   PetscPrintf(PETSC_COMM_WORLD,"vfix:%e \n",element->vfix[f]);
 //#############################################
*/

    }

    } else {

    for (f=0; f<n; f++) {
      PetscInt k = element->ifix[f];
      element->ufix[f] = U[k];
        //U[k] = element->vfix[f]*0.666666;
        V[k] = element->vfix[f];
        A[k] = 0.0;
    }

    }


  }
  PetscFunctionReturn(0);
}



#undef  __FUNCT__
#define __FUNCT__ "IGAElementFixValues"
PetscErrorCode IGAElementFixValues(IGAElement element,PetscScalar U[])
{
  PetscFunctionBegin;
  PetscValidPointer(element,1);
  PetscValidScalarPointer(U,2);
  {
    PetscInt f,n;
    n = element->nfix;
    for (f=0; f<n; f++) {
      PetscInt k = element->ifix[f];
      element->ufix[f] = U[k];
      U[k] = element->vfix[f];
/*
 //#############################################
   PetscPrintf(PETSC_COMM_WORLD,"ifix:%d \n",element->ifix[f]);
   PetscPrintf(PETSC_COMM_WORLD,"vfix:%e \n",element->vfix[f]);
 //#############################################
*/
    }
  }
  PetscFunctionReturn(0);
}



#undef  __FUNCT__
#define __FUNCT__ "IGAElementFixSystem"
PetscErrorCode IGAElementFixSystem(IGAElement element,PetscScalar K[],PetscScalar F[])
{
  PetscFunctionBegin;
  PetscValidPointer(element,1);
  PetscValidScalarPointer(K,2);
  PetscValidScalarPointer(F,3);
  if (PetscUnlikely(element->collocation)) goto collocation;
  {
    PetscInt M = element->neq * element->dof;
    PetscInt N = element->nen * element->dof;
    PetscInt f,n;
    n = element->nflux;
    for (f=0; f<n; f++) {
      PetscInt    k = element->iflux[f];
      PetscScalar v = element->vflux[f];
      F[k] += v;
    }
    n = element->nfix;
    for (f=0; f<n; f++) {
      PetscInt    k = element->ifix[f];
      PetscScalar v = element->vfix[f];
      PetscInt i,j;
      for (i=0; i<M; i++) F[i] -= K[i*N+k] * v;
      for (i=0; i<M; i++) K[i*N+k] = 0.0;
      for (j=0; j<N; j++) K[k*N+j] = 0.0;
      K[k*N+k] = 1.0;
      F[k]     = v;
    }
  }
  PetscFunctionReturn(0);
 collocation:
  {
    PetscInt dim = element->dim;
    PetscInt dof = element->dof;
    PetscInt nen = element->nen;
    PetscInt N = nen * dof;
    PetscInt dir,side;
    for (dir=0; dir<dim; dir++) {
      for (side=0; side<2; side++) {
        IGABoundary b = AtBoundary(element,dir,side);
        if(b && b->nload) {
          PetscInt  f, n = b->nload;
          PetscReal *dshape, normal[3] = {0.0,0.0,0.0};
          if (!element->geometry) {
            normal[dir] = side ? 1.0 : -1.0;
            dshape = element->basis[1];
          } else {
            PetscReal dS, *F = element->gradX[0];
            IGA_GetNormal(dim,dir,side,F,&dS,normal);
            dshape = element->shape[1];
          }
          for (f=0; f<n; f++) {
            PetscInt    c = b->iload[f];
            PetscScalar v = b->vload[f];
            PetscInt    a,j;
            for (j=0; j<N; j++) K[c*N+j] = 0.0;
            for (a=0; a<nen; a++)
              K[c*N+a*dof+c] = DOT(dim,&dshape[a*dim],normal);
            F[c] = v;
          }
        }
        if (b && b->count) {
          PetscInt  f, n = b->count;
          PetscReal *shape = element->basis[0];
          for (f=0; f<n; f++) {
            PetscInt    c = b->field[f];
            PetscScalar v = b->value[f];
            PetscInt    a,j;
            for (j=0; j<N; j++) K[c*N+j] = 0.0;
            for (a=0; a<nen; a++)
              K[c*N+a*dof+c] = shape[a];
            F[c] = v;
          }
        }
      }
    }
  }
  PetscFunctionReturn(0);
}



#undef  __FUNCT__
#define __FUNCT__ "IGAElementFixFunction_JB"
PetscErrorCode IGAElementFixFunction_JB(IGAElement element,PetscScalar F[])
{
  PetscFunctionBegin;
  PetscValidPointer(element,1);
  PetscValidScalarPointer(F,2);
  if (PetscUnlikely(element->collocation)) goto collocation;
  {
    PetscInt f,n;
    n = element->nflux;
    for (f=0; f<n; f++) {
      PetscInt    k = element->iflux[f];
      PetscScalar v = element->vflux[f];
      F[k] -= v;
    }
    n = element->nfix;
    for (f=0; f<n; f++) {
      PetscInt    k = element->ifix[f];
//      PetscScalar v = element->vfix[f];
//      PetscScalar u = element->ufix[f];
      F[k] = 0.0; //u-element->vfix[f]/(0.340277*1e-4);  // (u-element->vfix[f])*1.0/(0.44444444*1e-4);    //u-element->vfix[f]/(0.44444444*1e-4); //u - v;


/*
//#############################################
  PetscPrintf(PETSC_COMM_WORLD," ifix:%d ",element->ifix[f]);
  PetscPrintf(PETSC_COMM_WORLD," vfix:%e ",element->vfix[f]);
  PetscPrintf(PETSC_COMM_WORLD," ufix:%e ",element->ufix[f]);
  PetscPrintf(PETSC_COMM_WORLD," F(ifix):%e \n",F[k]);
//#############################################
*/

    }
  }
  PetscFunctionReturn(0);
 collocation:
  {
    PetscInt f,n;
    n = element->nfix;
    for (f=0; f<n; f++) {
      PetscInt k = element->ifix[f];
      PetscInt a = k / element->dof;
      PetscInt c = k % element->dof;
      if (element->rowmap[0] == element->colmap[a])
        {
          PetscScalar v = element->vfix[f];
          PetscScalar u = element->ufix[f];
          F[c] = u - v;
        }
    }
  }
  PetscFunctionReturn(0);
}






#undef  __FUNCT__
#define __FUNCT__ "IGAElementFixFunction"
PetscErrorCode IGAElementFixFunction(IGAElement element,PetscScalar F[])
{
  PetscFunctionBegin;
  PetscValidPointer(element,1);
  PetscValidScalarPointer(F,2);
  if (PetscUnlikely(element->collocation)) goto collocation;
  {
    PetscInt f,n;
    n = element->nflux;
    for (f=0; f<n; f++) {
      PetscInt    k = element->iflux[f];
      PetscScalar v = element->vflux[f];
      F[k] -= v;
    }
    n = element->nfix;
    for (f=0; f<n; f++) {
      PetscInt    k = element->ifix[f];
      PetscScalar v = element->vfix[f];
      PetscScalar u = element->ufix[f];
      F[k] = u - v;
/*
//#############################################
  PetscPrintf(PETSC_COMM_WORLD," ifix:%d ",element->ifix[f]);
  PetscPrintf(PETSC_COMM_WORLD," vfix:%e ",element->vfix[f]);
  PetscPrintf(PETSC_COMM_WORLD," ufix:%e ",element->ufix[f]);
  PetscPrintf(PETSC_COMM_WORLD," F(ifix):%e \n",F[k]);
//#############################################
*/
    }
  }
  PetscFunctionReturn(0);
 collocation:
  {
    PetscInt f,n;
    n = element->nfix;
    for (f=0; f<n; f++) {
      PetscInt k = element->ifix[f];
      PetscInt a = k / element->dof;
      PetscInt c = k % element->dof;
      if (element->rowmap[0] == element->colmap[a])
        {
          PetscScalar v = element->vfix[f];
          PetscScalar u = element->ufix[f];
          F[c] = u - v;
        }
    }
  }
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAElementFixJacobian"
PetscErrorCode IGAElementFixJacobian(IGAElement element,PetscScalar J[])
{
  PetscFunctionBegin;
  PetscValidPointer(element,1);
  PetscValidScalarPointer(J,2);
  if (PetscUnlikely(element->collocation)) goto collocation;
  {
    PetscInt M = element->neq * element->dof;
    PetscInt N = element->nen * element->dof;
    PetscInt f,n;
    n = element->nfix;
//    PetscPrintf(PETSC_COMM_WORLD," nen:%d",element->index);
//    PetscPrintf(PETSC_COMM_WORLD," nen:%d",element->nen);

    for (f=0; f<n; f++) {
      PetscInt i,j,k=element->ifix[f];
      for (i=0; i<M; i++) J[i*N+k] = 0.0;
      for (j=0; j<N; j++) J[k*N+j] = 0.0;
      J[k*N+k] = 1;

//      PetscPrintf(PETSC_COMM_WORLD," J(K*N+k):%d , %e \n",k, J[k*N+k]);

    }
  }
  PetscFunctionReturn(0);
 collocation:
  {
    PetscInt nen = element->nen;
    PetscInt dof = element->dof;
    PetscInt f,n;
    n = element->nfix;
    for (f=0; f<n; f++) {
      PetscInt k = element->ifix[f];
      PetscInt a = k / dof;
      PetscInt c = k % dof;
      if (element->rowmap[0] == element->colmap[a])
        {
          PetscInt  i,j,N=nen*dof;
          PetscReal *shape = element->basis[0];
          for (j=0; j<N; j++) J[c*N+j] = 0.0;
          for (i=0; i<nen; i++)
            J[c*N+i*dof+c] = shape[i];
        }
    }
  }
  PetscFunctionReturn(0);
}


#undef  __FUNCT__
#define __FUNCT__ "IGAElementFixJacobian_JB"
PetscErrorCode IGAElementFixJacobian_JB(IGAElement element,PetscScalar J[], TS ts)
{
  PetscFunctionBegin;
  PetscValidPointer(element,1);
  PetscValidScalarPointer(J,2);

  TS_Alpha2     *th    = (TS_Alpha2*)ts->data;

/*
  AppCtx *user = (AppCtx *)ctx;
  PetscPrintf(PETSC_COMM_WORLD,"dt %e \n", user->dt);
  PetscPrintf(PETSC_COMM_WORLD,"Beta %e \n", user->Beta);
  PetscPrintf(PETSC_COMM_WORLD,"Alpha_m %e \n", user->Alpha_m);
  PetscPrintf(PETSC_COMM_WORLD,"Alpha_f eeeeeeee %e \n", th->Alpha_f);
*/

  if (PetscUnlikely(element->collocation)) goto collocation;
  {
    PetscInt M = element->neq * element->dof;
    PetscInt N = element->nen * element->dof;
    PetscInt f,n;
    n = element->nfix;


    for (f=0; f<n; f++) {
      PetscInt i,j,k=element->ifix[f];
      for (i=0; i<M; i++) J[i*N+k] = 0.0;
      for (j=0; j<N; j++) J[k*N+j] = 0.0;
      J[k*N+k] = 1.0/(th->scale_J);
//      PetscPrintf(PETSC_COMM_WORLD," J(K*N+k):%d , %e \n",k, J[k*N+k]);

    }
  }
  PetscFunctionReturn(0);
 collocation:
  {
    PetscInt nen = element->nen;
    PetscInt dof = element->dof;
    PetscInt f,n;
    n = element->nfix;
    for (f=0; f<n; f++) {
      PetscInt k = element->ifix[f];
      PetscInt a = k / dof;
      PetscInt c = k % dof;
      if (element->rowmap[0] == element->colmap[a])
        {
          PetscInt  i,j,N=nen*dof;
          PetscReal *shape = element->basis[0];
          for (j=0; j<N; j++) J[c*N+j] = 0.0;
          for (i=0; i<nen; i++)
            J[c*N+i*dof+c] = shape[i];
        }
    }
  }
  PetscFunctionReturn(0);
}





#undef  __FUNCT__
#define __FUNCT__ "IGAElementAssembleVec"
PetscErrorCode IGAElementAssembleVec(IGAElement element,const PetscScalar F[],Vec vec)
{
  PetscInt       mm,*ii;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidPointer(element,1);
  PetscValidScalarPointer(F,2);
  PetscValidHeaderSpecific(vec,VEC_CLASSID,3);
  mm = element->neq; ii = element->rowmap;
  if (element->dof == 1) {
    ierr = VecSetValuesLocal(vec,mm,ii,F,ADD_VALUES);CHKERRQ(ierr);
  } else {
    ierr = VecSetValuesBlockedLocal(vec,mm,ii,F,ADD_VALUES);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAElementAssembleMat"
PetscErrorCode IGAElementAssembleMat(IGAElement element,const PetscScalar K[],Mat mat)
{
  PetscInt       mm,*ii;
  PetscInt       nn,*jj;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidPointer(element,1);
  PetscValidScalarPointer(K,2);
  PetscValidHeaderSpecific(mat,MAT_CLASSID,3);

  mm = element->neq; ii = element->rowmap;
  nn = element->nen; jj = element->colmap;
  if (element->dof == 1) {
    ierr = MatSetValuesLocal(mat,mm,ii,nn,jj,K,ADD_VALUES);CHKERRQ(ierr);
  } else {
    ierr = MatSetValuesBlockedLocal(mat,mm,ii,nn,jj,K,ADD_VALUES);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

/* -------------------------------------------------------------------------- */

EXTERN_C_BEGIN
extern void IGA_Interpolate(PetscInt nen,PetscInt dof,PetscInt dim,PetscInt der,
                            const PetscReal N[],const PetscScalar U[],PetscScalar u[]);
EXTERN_C_END

#define ELENGTH(a,b) sqrt((corners[a][0]-corners[b][0])*(corners[a][0]-corners[b][0])+(corners[a][1]-corners[b][1])*(corners[a][1]-corners[b][1])+(corners[a][2]-corners[b][2])*(corners[a][2]-corners[b][2]))

#undef  __FUNCT__
#define __FUNCT__ "IGAElementCharacteristicSize"
PetscErrorCode IGAElementCharacteristicSize(IGAElement element,PetscReal *h)
{
  PetscInt       dir,i,m,e,nc=0;
  PetscReal     *U;
  PetscReal      limits[3][2],corners[8][3],g[3]={0,0,0};
  IGAPoint       p;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidPointer(element,1);

  /* find parametric element breaks */
  for(dir=0;dir<element->dim;dir++){
    U = element->parent->axis[dir]->U;
    m = element->parent->axis[dir]->m;
    e = -1;
    for(i=0;i<m-1;i++){
      if(fabs(U[i+1]-U[i]) > 100*PETSC_MACHINE_EPSILON) e += 1;
      if(e==element->ID[dir]) { limits[dir][0] = U[i]; limits[dir][1] = U[i+1]; break;}
    }
  }

  /* simple subtraction if in parametric domain */
  *h = -1;
  if(!(element->geometry)){
    for(dir=0;dir<element->dim;dir++) *h = PetscMax(*h,limits[dir][1]-limits[dir][0]);
    PetscFunctionReturn(0);
  }

  /* find parametric corners */
  switch (element->dim) {
  case 3:
    nc = 8;
    corners[0][0] = limits[0][0]; corners[0][1] = limits[1][0]; corners[0][2] = limits[2][0];
    corners[1][0] = limits[0][1]; corners[1][1] = limits[1][0]; corners[1][2] = limits[2][0];
    corners[2][0] = limits[0][0]; corners[2][1] = limits[1][1]; corners[2][2] = limits[2][0];
    corners[3][0] = limits[0][1]; corners[3][1] = limits[1][1]; corners[3][2] = limits[2][0];
    corners[4][0] = limits[0][0]; corners[4][1] = limits[1][0]; corners[4][2] = limits[2][1];
    corners[5][0] = limits[0][1]; corners[5][1] = limits[1][0]; corners[5][2] = limits[2][1];
    corners[6][0] = limits[0][0]; corners[6][1] = limits[1][1]; corners[6][2] = limits[2][1];
    corners[7][0] = limits[0][1]; corners[7][1] = limits[1][1]; corners[7][2] = limits[2][1];
    break;
  case 2:
    nc = 4;
    corners[0][0] = limits[0][0]; corners[0][1] = limits[1][0]; corners[0][2] = 0;
    corners[1][0] = limits[0][1]; corners[1][1] = limits[1][0]; corners[1][2] = 0;
    corners[2][0] = limits[0][0]; corners[2][1] = limits[1][1]; corners[2][2] = 0;
    corners[3][0] = limits[0][1]; corners[3][1] = limits[1][1]; corners[3][2] = 0;
    break;
  case 1:
    nc = 2;
    corners[0][0] = limits[0][0]; corners[0][1] = 0; corners[0][2] = 0;
    corners[1][0] = limits[0][1]; corners[1][1] = 0; corners[1][2] = 0;
    break;
  }

  /* interpolate corners */
  ierr = IGAPointCreate(&p);CHKERRQ(ierr);
  for(i=0;i<nc;i++){
    p->point = &(corners[i][0]);
    ierr = IGAPointInit(p,element);CHKERRQ(ierr);
    ierr = IGAPointEval(element->parent,p);CHKERRQ(ierr);
    IGA_Interpolate(p->nen,p->dim,p->dim,0,p->basis[0],p->geometry,&g[0]);
    for(e=0;e<p->dim;e++) corners[i][e] = g[e];
  }
  ierr = IGAPointDestroy(&p);CHKERRQ(ierr);

  /* compute longest edge treating element as a quad/hex */
  switch (element->dim) {
  case 3:
    *h = PetscMax(*h,ELENGTH(1,0));
    *h = PetscMax(*h,ELENGTH(3,2));
    *h = PetscMax(*h,ELENGTH(0,2));
    *h = PetscMax(*h,ELENGTH(1,3));
    *h = PetscMax(*h,ELENGTH(4,5));
    *h = PetscMax(*h,ELENGTH(6,7));
    *h = PetscMax(*h,ELENGTH(4,6));
    *h = PetscMax(*h,ELENGTH(5,7));
    *h = PetscMax(*h,ELENGTH(0,4));
    *h = PetscMax(*h,ELENGTH(1,5));
    *h = PetscMax(*h,ELENGTH(2,6));
    *h = PetscMax(*h,ELENGTH(3,7));
    break;
  case 2:
    *h = PetscMax(*h,ELENGTH(1,0));
    *h = PetscMax(*h,ELENGTH(3,2));
    *h = PetscMax(*h,ELENGTH(0,2));
    *h = PetscMax(*h,ELENGTH(1,3));
    break;
  case 1:
    *h = PetscMax(*h,ELENGTH(1,0));
    break;
  }

  PetscFunctionReturn(0);
}
