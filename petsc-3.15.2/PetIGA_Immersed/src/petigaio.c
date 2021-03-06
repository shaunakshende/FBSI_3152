#include "petiga.h"
#include "petigagrid.h"

PETSC_EXTERN PetscErrorCode IGASetUp_Basic(IGA);
static       PetscErrorCode VecLoad_Binary_SkipHeader(Vec,PetscViewer);

#undef  __FUNCT__
#define __FUNCT__ "IGALoad"
PetscErrorCode IGALoad(IGA iga,PetscViewer viewer)
{
  PetscBool      isbinary;
  PetscBool      skipheader;
  PetscBool      geometry;
  PetscBool      property;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(iga,1,viewer,2);

  /* */
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERBINARY,&isbinary);CHKERRQ(ierr);
  if (!isbinary) SETERRQ(((PetscObject)viewer)->comm,PETSC_ERR_ARG_WRONG,"Only for binary viewers");
  ierr = PetscViewerBinaryGetSkipHeader(viewer,&skipheader);CHKERRQ(ierr);

  if (!skipheader) {
    PetscInt classid = 0;
    ierr = PetscViewerBinaryRead(viewer,&classid,1,PETSC_INT);CHKERRQ(ierr);
    if (classid != IGA_FILE_CLASSID) SETERRQ(((PetscObject)viewer)->comm,PETSC_ERR_ARG_WRONG,"Not an IGA in file");
  }
  { /* */
    PetscInt info = 0;
    ierr = PetscViewerBinaryRead(viewer,&info,1,PETSC_INT);CHKERRQ(ierr);
    geometry = (info & 0x1) ? PETSC_TRUE : PETSC_FALSE;
    property = (info & 0x2) ? PETSC_TRUE : PETSC_FALSE;
  }
  ierr = IGAReset(iga);CHKERRQ(ierr);
  { /* */
    PetscInt i,dim;
    ierr = PetscViewerBinaryRead(viewer,&dim, 1,PETSC_INT);CHKERRQ(ierr);
    ierr = IGASetDim(iga,dim);CHKERRQ(ierr);
    for (i=0; i<dim; i++) {
      IGAAxis   axis;
      PetscInt  p,m;
      PetscReal *U;
      ierr = PetscViewerBinaryRead(viewer,&p,1,PETSC_INT);CHKERRQ(ierr);
      ierr = PetscViewerBinaryRead(viewer,&m,1,PETSC_INT);CHKERRQ(ierr);
      ierr = PetscMalloc1(m,PetscReal,&U);CHKERRQ(ierr);
      ierr = PetscViewerBinaryRead(viewer,U,m,PETSC_REAL);CHKERRQ(ierr);
      ierr = IGAGetAxis(iga,i,&axis);CHKERRQ(ierr);
      ierr = IGAAxisInit(axis,p,m-1,U);CHKERRQ(ierr);CHKERRQ(ierr);
      ierr = PetscFree(U);CHKERRQ(ierr);
    }
  }
  ierr = IGASetUp_Basic(iga);CHKERRQ(ierr);
  if (geometry) { /* */
    PetscInt dim;
    ierr = PetscViewerBinaryRead(viewer,&dim,1,PETSC_INT);CHKERRQ(ierr);
    ierr = IGASetGeometryDim(iga,dim);CHKERRQ(ierr);
    ierr = IGALoadGeometry(iga,viewer);CHKERRQ(ierr);
  }
  if (property) { /* */
    PetscInt dim;
    ierr = PetscViewerBinaryRead(viewer,&dim,1,PETSC_INT);CHKERRQ(ierr);
    ierr = IGASetPropertyDim(iga,dim);CHKERRQ(ierr);
    ierr = IGALoadProperty(iga,viewer);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGASave"
PetscErrorCode IGASave(IGA iga,PetscViewer viewer)
{
  PetscBool      isbinary;
  PetscBool      skipheader;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  IGACheckSetUp(iga,1);

  if (viewer) {
    PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
    PetscCheckSameComm(iga,1,viewer,2);
  } else {
    MPI_Comm comm;
    ierr = IGAGetComm(iga,&comm);CHKERRQ(ierr);
    viewer = PETSC_VIEWER_BINARY_(comm);
    if (!viewer) PetscFunctionReturn(PETSC_ERR_PLIB);
  }
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERBINARY,&isbinary);CHKERRQ(ierr);
  if (!isbinary) SETERRQ(((PetscObject)viewer)->comm,PETSC_ERR_ARG_WRONG,"Only for binary viewers");
  ierr = PetscViewerBinaryGetSkipHeader(viewer,&skipheader);CHKERRQ(ierr);

  if (!skipheader) {
    PetscInt classid = IGA_FILE_CLASSID;
    ierr = PetscViewerBinaryWrite(viewer,&classid,1,PETSC_INT,PETSC_TRUE);CHKERRQ(ierr);
  }
  { /* */
    PetscInt info = 0;
    if (iga->geometry) info |= 0x1;
    if (iga->property) info |= 0x2;
    ierr = PetscViewerBinaryWrite(viewer,&info,1,PETSC_INT,PETSC_TRUE);CHKERRQ(ierr);
  }
  { /* */
    PetscInt i,dim;
    ierr = IGAGetDim(iga,&dim);CHKERRQ(ierr);
    ierr = PetscViewerBinaryWrite(viewer,&dim,1,PETSC_INT,PETSC_FALSE);CHKERRQ(ierr);
    for (i=0; i<dim; i++) {
      IGAAxis   axis;
      PetscInt  p,m,buf[2];
      PetscReal *U;
      ierr = IGAGetAxis(iga,i,&axis);CHKERRQ(ierr);
      ierr = IGAAxisGetDegree(axis,&p);CHKERRQ(ierr);
      ierr = IGAAxisGetKnots(axis,&m,&U);CHKERRQ(ierr);
      buf[0] = p; buf[1] = m+1;
      ierr = PetscViewerBinaryWrite(viewer,buf,2,PETSC_INT,PETSC_TRUE);CHKERRQ(ierr);
      ierr = PetscViewerBinaryWrite(viewer,U,m+1,PETSC_REAL,PETSC_FALSE);CHKERRQ(ierr);
    }
  }
  if (iga->geometry) { /* */
    PetscInt dim;
    ierr = IGAGetGeometryDim(iga,&dim);CHKERRQ(ierr);
    ierr = PetscViewerBinaryWrite(viewer,&dim,1,PETSC_INT,PETSC_TRUE);CHKERRQ(ierr);
    ierr = IGASaveGeometry(iga,viewer);CHKERRQ(ierr);
  }
  if (iga->property) { /* */
    PetscInt dim;
    ierr = IGAGetPropertyDim(iga,&dim);CHKERRQ(ierr);
    ierr = PetscViewerBinaryWrite(viewer,&dim,1,PETSC_INT,PETSC_TRUE);CHKERRQ(ierr);
    ierr = IGASaveProperty(iga,viewer);CHKERRQ(ierr);
  }
  ierr = PetscViewerFlush(viewer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

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

#undef  __FUNCT__
#define __FUNCT__ "IGASetGeometryDim"
/*@
   IGASetGeometryDim - Sets the dimension of the geometry

   Logically Collective on IGA

   Input Parameters:
+  iga - the IGA context
-  dim - the dimension of the geometry

   Level: normal

.keywords: IGA, dimension
@*/
PetscErrorCode IGASetGeometryDim(IGA iga,PetscInt dim)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidLogicalCollectiveInt(iga,dim,2);
  if (dim < 1 || dim > 3)
    SETERRQ1(((PetscObject)iga)->comm,PETSC_ERR_ARG_WRONGSTATE,
             "Number of space dimensions must be in range [1,3], got %D",dim);
  if (iga->geometry > 0 && iga->geometry != dim)
    SETERRQ2(((PetscObject)iga)->comm,PETSC_ERR_ARG_WRONGSTATE,
             "Cannot change IGA spatial dim to %D after it was set to %D",dim,iga->geometry);
  if (iga->geometry == 0) iga->setup = PETSC_FALSE;
  iga->geometry = dim;
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAGetGeometryDim"
PetscErrorCode IGAGetGeometryDim(IGA iga,PetscInt *dim)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidPointer(dim,2);
  *dim = iga->geometry;
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGALoadGeometry"
PetscErrorCode IGALoadGeometry(IGA iga,PetscViewer viewer)
{
  PetscBool      isbinary;
  PetscBool      skipheader;
  PetscInt       nsd;
  PetscReal      min_w,max_w,tol_w = 100*PETSC_MACHINE_EPSILON;
  Vec            nvec,gvec,lvec;
  VecScatter     g2n,g2l;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(iga,1,viewer,2);
  if (iga->setupstage < 1) IGACheckSetUp(iga,1);

  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERBINARY,&isbinary);CHKERRQ(ierr);
  if (!isbinary) SETERRQ(((PetscObject)viewer)->comm,PETSC_ERR_ARG_WRONG,"Only for binary viewers");
  ierr = PetscViewerBinaryGetSkipHeader(viewer,&skipheader);CHKERRQ(ierr);

  ierr = IGAGetGeometryDim(iga,&nsd);CHKERRQ(ierr);
  if (nsd < 1) SETERRQ(((PetscObject)iga)->comm,PETSC_ERR_ARG_WRONGSTATE,
                       "Must call IGASetGeometryDim() first");
  {
    IGA_Grid grid;
    ierr = IGA_NewGridIO(iga,nsd+1,&grid);CHKERRQ(ierr);
    ierr = IGA_Grid_GetVecNatural(grid,iga->vectype,&nvec);CHKERRQ(ierr);
    ierr = IGA_Grid_GetVecGlobal (grid,iga->vectype,&gvec);CHKERRQ(ierr);
    ierr = IGA_Grid_GetVecLocal  (grid,iga->vectype,&lvec);CHKERRQ(ierr);
    ierr = IGA_Grid_GetScatterG2N(grid,&g2n);CHKERRQ(ierr);
    ierr = IGA_Grid_GetScatterG2L(grid,&g2l);CHKERRQ(ierr);
    ierr = PetscObjectReference((PetscObject)nvec);CHKERRQ(ierr);
    ierr = PetscObjectReference((PetscObject)gvec);CHKERRQ(ierr);
    ierr = PetscObjectReference((PetscObject)lvec);CHKERRQ(ierr);
    ierr = PetscObjectReference((PetscObject)g2n);CHKERRQ(ierr);
    ierr = PetscObjectReference((PetscObject)g2l);CHKERRQ(ierr);
    ierr = IGA_Grid_Destroy(&grid);CHKERRQ(ierr);
  }
  /* viewer -> natural*/
  if (!skipheader)
    {ierr = VecLoad(nvec,viewer);CHKERRQ(ierr);}
  else
    {ierr = VecLoad_Binary_SkipHeader(nvec,viewer);CHKERRQ(ierr);}
  /* natural -> global */
  ierr = VecScatterBegin(g2n,nvec,gvec,INSERT_VALUES,SCATTER_REVERSE);CHKERRQ(ierr);
  ierr = VecScatterEnd  (g2n,nvec,gvec,INSERT_VALUES,SCATTER_REVERSE);CHKERRQ(ierr);
  /* global -> local */
  ierr = VecScatterBegin(g2l,gvec,lvec,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  ierr = VecScatterEnd  (g2l,gvec,lvec,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  ierr = VecStrideMin(gvec,nsd,NULL,&min_w);CHKERRQ(ierr);
  ierr = VecStrideMax(gvec,nsd,NULL,&max_w);CHKERRQ(ierr);

  iga->rational = ((max_w-min_w)>tol_w) ? PETSC_TRUE : PETSC_FALSE;
  ierr = PetscFree(iga->geometryX);CHKERRQ(ierr);
  ierr = PetscFree(iga->geometryINI);CHKERRQ(ierr);
  ierr = PetscFree(iga->rationalW);CHKERRQ(ierr);
   {
    PetscInt n,a,i,pos;
    PetscReal *X,*W, *INI;
    const PetscScalar *Xw;
    ierr = VecGetSize(lvec,&n);CHKERRQ(ierr);
    n /= (nsd+1);
    ierr = PetscMalloc1(n*nsd,PetscReal,&iga->geometryX);CHKERRQ(ierr);
    ierr = PetscMalloc1(n*nsd,PetscReal,&iga->geometryINI);CHKERRQ(ierr);
    ierr = PetscMalloc1(n,    PetscReal,&iga->rationalW);CHKERRQ(ierr);
    X = iga->geometryX; W = iga->rationalW; INI = iga->geometryINI;
    ierr = VecGetArrayRead(lvec,&Xw);CHKERRQ(ierr);
    for (pos=0,a=0; a<n; a++) {
      for (i=0; i<nsd; i++) {
       X[i+a*nsd]   = PetscRealPart(Xw[pos++]);
       INI[i+a*nsd] = X[i+a*nsd];
      }
      W[a] = PetscRealPart(Xw[pos++]);
      if (W[a] != 0.0)
        for (i=0; i<nsd; i++) {
          X[i+a*nsd]   /= W[a];
          INI[i+a*nsd] /= W[a];
        }
    }
    ierr = VecRestoreArrayRead(lvec,&Xw);CHKERRQ(ierr);
  }

  ierr = VecScatterDestroy(&g2n);CHKERRQ(ierr);
  ierr = VecScatterDestroy(&g2l);CHKERRQ(ierr);
  ierr = VecDestroy(&lvec);CHKERRQ(ierr);
  ierr = VecDestroy(&gvec);CHKERRQ(ierr);
  ierr = VecDestroy(&nvec);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGASaveGeometry"
PetscErrorCode IGASaveGeometry(IGA iga,PetscViewer viewer)
{
  PetscBool      isbinary;
  PetscBool      skipheader;
  PetscInt       nsd;
  Vec            nvec,gvec,lvec;
  VecScatter     l2g,g2n;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(iga,1,viewer,2);
  if (iga->setupstage < 1) IGACheckSetUp(iga,1);
  if (!iga->geometry) SETERRQ(((PetscObject)iga)->comm,PETSC_ERR_ARG_WRONGSTATE,"No geometry set");

  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERBINARY,&isbinary);CHKERRQ(ierr);
  if (!isbinary) SETERRQ(((PetscObject)viewer)->comm,PETSC_ERR_ARG_WRONG,"Only for binary viewers");
  ierr = PetscViewerBinaryGetSkipHeader(viewer,&skipheader);CHKERRQ(ierr);

  ierr = IGAGetGeometryDim(iga,&nsd);CHKERRQ(ierr);
  {
    IGA_Grid grid;
    ierr = IGA_NewGridIO(iga,nsd+1,&grid);CHKERRQ(ierr);
    ierr = IGA_Grid_GetVecNatural(grid,iga->vectype,&nvec);CHKERRQ(ierr);
    ierr = IGA_Grid_GetVecGlobal (grid,iga->vectype,&gvec);CHKERRQ(ierr);
    ierr = IGA_Grid_GetVecLocal  (grid,iga->vectype,&lvec);CHKERRQ(ierr);
    ierr = IGA_Grid_GetScatterL2G(grid,&l2g);CHKERRQ(ierr);
    ierr = IGA_Grid_GetScatterG2N(grid,&g2n);CHKERRQ(ierr);
    ierr = PetscObjectReference((PetscObject)nvec);CHKERRQ(ierr);
    ierr = PetscObjectReference((PetscObject)gvec);CHKERRQ(ierr);
    ierr = PetscObjectReference((PetscObject)lvec);CHKERRQ(ierr);
    ierr = PetscObjectReference((PetscObject)l2g);CHKERRQ(ierr);
    ierr = PetscObjectReference((PetscObject)g2n);CHKERRQ(ierr);
    ierr = IGA_Grid_Destroy(&grid);CHKERRQ(ierr);
  }
  {
    PetscInt n,a,i,pos;
    PetscScalar *Xw;
    const PetscReal *X = iga->geometryX;
    const PetscReal *W = iga->rationalW;
    ierr = VecGetSize(lvec,&n);CHKERRQ(ierr);
    n /= (nsd+1);
    ierr = VecGetArray(lvec,&Xw);CHKERRQ(ierr);
    for (pos=0,a=0; a<n; a++) {
      PetscReal w = (W[a] != 0.0) ? W[a] : 1.0;
      for (i=0; i<nsd; i++)
        Xw[pos++] = X[i+a*nsd] * w;
      Xw[pos++] = W[a];
    }
    ierr = VecRestoreArray(lvec,&Xw);CHKERRQ(ierr);
  }
  /* local -> global */
  ierr = VecScatterBegin(l2g,lvec,gvec,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  ierr = VecScatterEnd  (l2g,lvec,gvec,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  /* global -> natural */
  ierr = VecScatterBegin(g2n,gvec,nvec,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  ierr = VecScatterEnd  (g2n,gvec,nvec,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  /* natural -> viewer */
  ierr = VecView(nvec,viewer);CHKERRQ(ierr);

  ierr = VecScatterDestroy(&g2n);CHKERRQ(ierr);
  ierr = VecScatterDestroy(&l2g);CHKERRQ(ierr);
  ierr = VecDestroy(&lvec);CHKERRQ(ierr);
  ierr = VecDestroy(&gvec);CHKERRQ(ierr);
  ierr = VecDestroy(&nvec);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGASetPropertyDim"
/*@
   IGASetPropertyDim - Sets the dimension of the property

   Logically Collective on IGA

   Input Parameters:
+  iga - the IGA context
-  dim - the dimension of the property

   Level: normal

.keywords: IGA, dimension
@*/
PetscErrorCode IGASetPropertyDim(IGA iga,PetscInt dim)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidLogicalCollectiveInt(iga,dim,2);
  if (dim < 1)
    SETERRQ1(((PetscObject)iga)->comm,PETSC_ERR_ARG_WRONGSTATE,
             "Number of properties must be greater than 0, got %D",dim);
  if (iga->property > 0 && iga->property != dim)
    SETERRQ2(((PetscObject)iga)->comm,PETSC_ERR_ARG_WRONGSTATE,
             "Cannot change IGA property dim to %D after it was set to %D",dim,iga->property);
  if (iga->property == 0) iga->setup = PETSC_FALSE;
  iga->property = dim;
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAGetPropertyDim"
PetscErrorCode IGAGetPropertyDim(IGA iga,PetscInt *dim)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidPointer(dim,2);
  *dim = iga->property;
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGALoadProperty"
PetscErrorCode IGALoadProperty(IGA iga,PetscViewer viewer)
{
  PetscBool      isbinary;
  PetscBool      skipheader;
  PetscInt       npd;
  Vec            nvec,gvec,lvec;
  VecScatter     g2n,g2l;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(iga,1,viewer,2);
  if (iga->setupstage < 1) IGACheckSetUp(iga,1);

  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERBINARY,&isbinary);CHKERRQ(ierr);
  if (!isbinary) SETERRQ(((PetscObject)viewer)->comm,PETSC_ERR_ARG_WRONG,"Only for binary viewers");
  ierr = PetscViewerBinaryGetSkipHeader(viewer,&skipheader);CHKERRQ(ierr);

  ierr = IGAGetPropertyDim(iga,&npd);CHKERRQ(ierr);
  if (npd < 1) SETERRQ(((PetscObject)iga)->comm,PETSC_ERR_ARG_WRONGSTATE,
                       "Must call IGASetPropertyDim() first");
  {
    IGA_Grid grid;
    ierr = IGA_NewGridIO(iga,npd,&grid);CHKERRQ(ierr);
    ierr = IGA_Grid_GetVecNatural(grid,iga->vectype,&nvec);CHKERRQ(ierr);
    ierr = IGA_Grid_GetVecGlobal (grid,iga->vectype,&gvec);CHKERRQ(ierr);
    ierr = IGA_Grid_GetVecLocal  (grid,iga->vectype,&lvec);CHKERRQ(ierr);
    ierr = IGA_Grid_GetScatterG2N(grid,&g2n);CHKERRQ(ierr);
    ierr = IGA_Grid_GetScatterG2L(grid,&g2l);CHKERRQ(ierr);
    ierr = PetscObjectReference((PetscObject)nvec);CHKERRQ(ierr);
    ierr = PetscObjectReference((PetscObject)gvec);CHKERRQ(ierr);
    ierr = PetscObjectReference((PetscObject)lvec);CHKERRQ(ierr);
    ierr = PetscObjectReference((PetscObject)g2n);CHKERRQ(ierr);
    ierr = PetscObjectReference((PetscObject)g2l);CHKERRQ(ierr);
    ierr = IGA_Grid_Destroy(&grid);CHKERRQ(ierr);
  }
  /* viewer -> natural*/
  if (!skipheader)
    {ierr = VecLoad(nvec,viewer);CHKERRQ(ierr);}
  else
    {ierr = VecLoad_Binary_SkipHeader(nvec,viewer);CHKERRQ(ierr);}
  /* natural -> global */
  ierr = VecScatterBegin(g2n,nvec,gvec,INSERT_VALUES,SCATTER_REVERSE);CHKERRQ(ierr);
  ierr = VecScatterEnd  (g2n,nvec,gvec,INSERT_VALUES,SCATTER_REVERSE);CHKERRQ(ierr);
  /* global -> local */
  ierr = VecScatterBegin(g2l,gvec,lvec,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  ierr = VecScatterEnd  (g2l,gvec,lvec,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);

  ierr = PetscFree(iga->propertyA);CHKERRQ(ierr);
  {
    PetscInt n; const PetscScalar *A;
    ierr = VecGetSize(lvec,&n);CHKERRQ(ierr);
    ierr = PetscMalloc1(n,PetscScalar,&iga->propertyA);CHKERRQ(ierr);
    ierr = VecGetArrayRead(lvec,&A);CHKERRQ(ierr);
    ierr = PetscMemcpy(iga->propertyA,A,n*sizeof(PetscScalar));CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(lvec,&A);CHKERRQ(ierr);
  }

  ierr = VecScatterDestroy(&g2n);CHKERRQ(ierr);
  ierr = VecScatterDestroy(&g2l);CHKERRQ(ierr);
  ierr = VecDestroy(&lvec);CHKERRQ(ierr);
  ierr = VecDestroy(&gvec);CHKERRQ(ierr);
  ierr = VecDestroy(&nvec);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGASaveProperty"
PetscErrorCode IGASaveProperty(IGA iga,PetscViewer viewer)
{
  PetscBool      isbinary;
  PetscBool      skipheader;
  PetscInt       npd;
  Vec            nvec,gvec,lvec;
  VecScatter     l2g,g2n;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(iga,1,viewer,2);
  if (iga->setupstage < 1) IGACheckSetUp(iga,1);
  if (!iga->property) SETERRQ(((PetscObject)iga)->comm,PETSC_ERR_ARG_WRONGSTATE,"No property set");

  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERBINARY,&isbinary);CHKERRQ(ierr);
  if (!isbinary) SETERRQ(((PetscObject)viewer)->comm,PETSC_ERR_ARG_WRONG,"Only for binary viewers");
  ierr = PetscViewerBinaryGetSkipHeader(viewer,&skipheader);CHKERRQ(ierr);

  ierr = IGAGetPropertyDim(iga,&npd);CHKERRQ(ierr);
  {
    IGA_Grid grid;
    ierr = IGA_NewGridIO(iga,npd,&grid);CHKERRQ(ierr);
    ierr = IGA_Grid_GetVecNatural(grid,iga->vectype,&nvec);CHKERRQ(ierr);
    ierr = IGA_Grid_GetVecGlobal (grid,iga->vectype,&gvec);CHKERRQ(ierr);
    ierr = IGA_Grid_GetVecLocal  (grid,iga->vectype,&lvec);CHKERRQ(ierr);
    ierr = IGA_Grid_GetScatterL2G(grid,&l2g);CHKERRQ(ierr);
    ierr = IGA_Grid_GetScatterG2N(grid,&g2n);CHKERRQ(ierr);
    ierr = PetscObjectReference((PetscObject)nvec);CHKERRQ(ierr);
    ierr = PetscObjectReference((PetscObject)gvec);CHKERRQ(ierr);
    ierr = PetscObjectReference((PetscObject)lvec);CHKERRQ(ierr);
    ierr = PetscObjectReference((PetscObject)l2g);CHKERRQ(ierr);
    ierr = PetscObjectReference((PetscObject)g2n);CHKERRQ(ierr);
    ierr = IGA_Grid_Destroy(&grid);CHKERRQ(ierr);
  }
  {
    PetscInt n; PetscScalar *A;
    ierr = VecGetSize(lvec,&n);CHKERRQ(ierr);
    ierr = VecGetArray(lvec,&A);CHKERRQ(ierr);
    ierr = PetscMemcpy(A,iga->propertyA,n*sizeof(PetscScalar));CHKERRQ(ierr);
    ierr = VecRestoreArray(lvec,&A);CHKERRQ(ierr);
  }
  /* local -> global */
  ierr = VecScatterBegin(l2g,lvec,gvec,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  ierr = VecScatterEnd  (l2g,lvec,gvec,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  /* global -> natural */
  ierr = VecScatterBegin(g2n,gvec,nvec,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  ierr = VecScatterEnd  (g2n,gvec,nvec,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  /* natural -> viewer */
  ierr = VecView(nvec,viewer);CHKERRQ(ierr);

  ierr = VecScatterDestroy(&g2n);CHKERRQ(ierr);
  ierr = VecScatterDestroy(&l2g);CHKERRQ(ierr);
  ierr = VecDestroy(&lvec);CHKERRQ(ierr);
  ierr = VecDestroy(&gvec);CHKERRQ(ierr);
  ierr = VecDestroy(&nvec);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGARead"
/*@
   IGARead - reads a IGA which has been saved in binary format

   Collective on IGA

   Input Parameters:
+  iga - the IGA context
-  filename - the file name which contains the IGA information

   Level: normal

.keywords: IGA, read
@*/
PetscErrorCode IGARead(IGA iga,const char filename[])
{
  MPI_Comm       comm;
  PetscViewer    viewer;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidCharPointer(filename,2);
  ierr = IGAGetComm(iga,&comm);CHKERRQ(ierr);
  ierr = PetscViewerCreate(comm,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer,PETSCVIEWERBINARY);CHKERRQ(ierr);
  ierr = PetscViewerBinarySkipInfo(viewer);CHKERRQ(ierr);
  /*ierr = PetscViewerBinarySetSkipHeader(viewer,PETSC_TRUE);CHKERRQ(ierr);*/
  ierr = PetscViewerFileSetMode(viewer,FILE_MODE_READ);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer,filename);CHKERRQ(ierr);
  ierr = IGALoad(iga,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAWrite"
/*@
   IGAWrite - writes a IGA to a file in binary format

   Collective on IGA

   Input Parameters:
+  iga - the IGA context
-  filename - the file name in which the IGA information is saved

   Level: normal

.keywords: IGA, write
@*/
PetscErrorCode IGAWrite(IGA iga,const char filename[])
{
  MPI_Comm       comm;
  PetscViewer    viewer;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidCharPointer(filename,2);
  ierr = IGAGetComm(iga,&comm);CHKERRQ(ierr);
  ierr = PetscViewerCreate(comm,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer,PETSCVIEWERBINARY);CHKERRQ(ierr);
  ierr = PetscViewerBinarySkipInfo(viewer);CHKERRQ(ierr);
  /*ierr = PetscViewerBinarySetSkipHeader(viewer,PETSC_TRUE);CHKERRQ(ierr);*/
  ierr = PetscViewerFileSetMode(viewer,FILE_MODE_WRITE);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer,filename);CHKERRQ(ierr);
  ierr = IGASave(iga,viewer);CHKERRQ(ierr);
  ierr = PetscViewerFlush(viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef  __FUNCT__
#define __FUNCT__ "VecLoad_Binary_SkipHeader"
static PetscErrorCode VecLoad_Binary_SkipHeader(Vec vec, PetscViewer viewer)
{
  MPI_Comm       comm;
  PetscMPIInt    i,rank,size,tag;
  int            fd;
  PetscInt       n;
  const PetscInt *range;
  PetscScalar    *array,*work;
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr = PetscObjectGetComm((PetscObject)viewer,&comm);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);
  ierr = PetscCommGetNewTag(comm,&tag);CHKERRQ(ierr);
  ierr = PetscViewerBinaryGetDescriptor(viewer,&fd);CHKERRQ(ierr);

  ierr = VecGetLocalSize(vec,&n);CHKERRQ(ierr);
  ierr = VecGetArray(vec,&array);CHKERRQ(ierr);
  if (!rank) {
    ierr = PetscBinaryRead(fd,array,n,PETSC_SCALAR);CHKERRQ(ierr);
    if (size > 1) {
      ierr = VecGetOwnershipRanges(vec,&range);CHKERRQ(ierr);
      n = 1;
      for (i=1; i<size; i++)
        n = PetscMax(n,range[i+1] - range[i]);
      ierr = PetscMalloc(n*sizeof(PetscScalar),&work);CHKERRQ(ierr);
      for (i=1; i<size; i++) {
        n   = range[i+1] - range[i];
        ierr = PetscBinaryRead(fd,work,n,PETSC_SCALAR);CHKERRQ(ierr);
        ierr = MPI_Send(work,n,MPIU_SCALAR,i,tag,comm);CHKERRQ(ierr);
      }
      ierr = PetscFree(work);CHKERRQ(ierr);
    }
  } else {
    MPI_Status status;
    ierr = MPI_Recv(array,n,MPIU_SCALAR,0,tag,comm,&status);CHKERRQ(ierr);
  }
  ierr = VecRestoreArray(vec,&array);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGALoadVec"
PetscErrorCode IGALoadVec(IGA iga,Vec vec,PetscViewer viewer)
{
  Vec            natural;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidHeaderSpecific(vec,VEC_CLASSID,2);
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,3);
  PetscCheckSameComm(iga,1,vec,2);
  PetscCheckSameComm(iga,1,viewer,3);
  IGACheckSetUp(iga,1);

  ierr = IGAGetNaturalVec(iga,&natural);
  ierr = VecLoad(natural,viewer);CHKERRQ(ierr);
  ierr = IGANaturalToGlobal(iga,natural,vec);

  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGASaveVec"
PetscErrorCode IGASaveVec(IGA iga,Vec vec,PetscViewer viewer)
{
  Vec            natural;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidHeaderSpecific(vec,VEC_CLASSID,2);
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,3);
  PetscCheckSameComm(iga,1,vec,2);
  PetscCheckSameComm(iga,1,viewer,3);
  IGACheckSetUp(iga,1);

  ierr = IGAGetNaturalVec(iga,&natural);
  ierr = IGAGlobalToNatural(iga,vec,natural);
  ierr = PetscObjectSetName((PetscObject)natural,((PetscObject)vec)->name);CHKERRQ(ierr);
  ierr = VecView(natural,viewer);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAReadVec"
PetscErrorCode IGAReadVec(IGA iga,Vec vec,const char filename[])
{
  MPI_Comm       comm;
  PetscViewer    viewer;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidHeaderSpecific(vec,VEC_CLASSID,1);
  PetscCheckSameComm(iga,1,vec,2);
  PetscValidCharPointer(filename,2);
  IGACheckSetUp(iga,1);

  ierr = IGAGetComm(iga,&comm);CHKERRQ(ierr);
  ierr = PetscViewerCreate(comm,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer,PETSCVIEWERBINARY);CHKERRQ(ierr);
  ierr = PetscViewerBinarySkipInfo(viewer);CHKERRQ(ierr);
  /*ierr = PetscViewerBinarySetSkipHeader(viewer,PETSC_TRUE);CHKERRQ(ierr);*/
  ierr = PetscViewerFileSetMode(viewer,FILE_MODE_READ);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer,filename);CHKERRQ(ierr);
  ierr = IGALoadVec(iga,vec,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAWriteVec"
PetscErrorCode IGAWriteVec(IGA iga,Vec vec,const char filename[])
{
  MPI_Comm       comm;
  PetscViewer    viewer;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidHeaderSpecific(vec,VEC_CLASSID,1);
  PetscCheckSameComm(iga,1,vec,2);
  PetscValidCharPointer(filename,2);
  IGACheckSetUp(iga,1);

  ierr = IGAGetComm(iga,&comm);CHKERRQ(ierr);
  ierr = PetscViewerCreate(comm,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer,PETSCVIEWERBINARY);CHKERRQ(ierr);
  ierr = PetscViewerBinarySkipInfo(viewer);CHKERRQ(ierr);
  /*ierr = PetscViewerBinarySetSkipHeader(viewer,PETSC_TRUE);CHKERRQ(ierr);*/
  ierr = PetscViewerFileSetMode(viewer,FILE_MODE_WRITE);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer,filename);CHKERRQ(ierr);
  ierr = IGASaveVec(iga,vec,viewer);CHKERRQ(ierr);
  ierr = PetscViewerFlush(viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
