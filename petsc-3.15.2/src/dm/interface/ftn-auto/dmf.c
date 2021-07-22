#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* dm.c */
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

#include "petscdm.h"
#include "petscdmlabel.h"
#include "petscds.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmcreate_ DMCREATE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmcreate_ dmcreate
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmclone_ DMCLONE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmclone_ dmclone
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define vecgetdm_ VECGETDM
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define vecgetdm_ vecgetdm
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define vecsetdm_ VECSETDM
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define vecsetdm_ vecsetdm
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define matgetdm_ MATGETDM
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define matgetdm_ matgetdm
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define matsetdm_ MATSETDM
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define matsetdm_ matsetdm
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmsetup_ DMSETUP
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmsetup_ dmsetup
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmsetfromoptions_ DMSETFROMOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmsetfromoptions_ dmsetfromoptions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmcreateglobalvector_ DMCREATEGLOBALVECTOR
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmcreateglobalvector_ dmcreateglobalvector
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmcreatelocalvector_ DMCREATELOCALVECTOR
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmcreatelocalvector_ dmcreatelocalvector
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetlocaltoglobalmapping_ DMGETLOCALTOGLOBALMAPPING
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetlocaltoglobalmapping_ dmgetlocaltoglobalmapping
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetblocksize_ DMGETBLOCKSIZE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetblocksize_ dmgetblocksize
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmcreateinterpolationscale_ DMCREATEINTERPOLATIONSCALE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmcreateinterpolationscale_ dmcreateinterpolationscale
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmcreaterestriction_ DMCREATERESTRICTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmcreaterestriction_ dmcreaterestriction
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmcreateinjection_ DMCREATEINJECTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmcreateinjection_ dmcreateinjection
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmcreatemassmatrix_ DMCREATEMASSMATRIX
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmcreatemassmatrix_ dmcreatemassmatrix
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmcreatecoloring_ DMCREATECOLORING
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmcreatecoloring_ dmcreatecoloring
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmcreatematrix_ DMCREATEMATRIX
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmcreatematrix_ dmcreatematrix
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmsetmatrixpreallocateonly_ DMSETMATRIXPREALLOCATEONLY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmsetmatrixpreallocateonly_ dmsetmatrixpreallocateonly
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmsetmatrixstructureonly_ DMSETMATRIXSTRUCTUREONLY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmsetmatrixstructureonly_ dmsetmatrixstructureonly
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmcreatesubdm_ DMCREATESUBDM
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmcreatesubdm_ dmcreatesubdm
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmrefine_ DMREFINE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmrefine_ dmrefine
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dminterpolate_ DMINTERPOLATE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dminterpolate_ dminterpolate
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dminterpolatesolution_ DMINTERPOLATESOLUTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dminterpolatesolution_ dminterpolatesolution
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetrefinelevel_ DMGETREFINELEVEL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetrefinelevel_ dmgetrefinelevel
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmsetrefinelevel_ DMSETREFINELEVEL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmsetrefinelevel_ dmsetrefinelevel
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmhasbasistransform_ DMHASBASISTRANSFORM
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmhasbasistransform_ dmhasbasistransform
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmglobaltolocal_ DMGLOBALTOLOCAL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmglobaltolocal_ dmglobaltolocal
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmglobaltolocalbegin_ DMGLOBALTOLOCALBEGIN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmglobaltolocalbegin_ dmglobaltolocalbegin
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmglobaltolocalend_ DMGLOBALTOLOCALEND
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmglobaltolocalend_ dmglobaltolocalend
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmlocaltoglobal_ DMLOCALTOGLOBAL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmlocaltoglobal_ dmlocaltoglobal
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmlocaltoglobalbegin_ DMLOCALTOGLOBALBEGIN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmlocaltoglobalbegin_ dmlocaltoglobalbegin
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmlocaltoglobalend_ DMLOCALTOGLOBALEND
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmlocaltoglobalend_ dmlocaltoglobalend
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmlocaltolocalbegin_ DMLOCALTOLOCALBEGIN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmlocaltolocalbegin_ dmlocaltolocalbegin
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmlocaltolocalend_ DMLOCALTOLOCALEND
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmlocaltolocalend_ dmlocaltolocalend
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmcoarsen_ DMCOARSEN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmcoarsen_ dmcoarsen
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmrestrict_ DMRESTRICT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmrestrict_ dmrestrict
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmsubdomainrestrict_ DMSUBDOMAINRESTRICT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmsubdomainrestrict_ dmsubdomainrestrict
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetcoarsenlevel_ DMGETCOARSENLEVEL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetcoarsenlevel_ dmgetcoarsenlevel
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmsetcoarsenlevel_ DMSETCOARSENLEVEL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmsetcoarsenlevel_ dmsetcoarsenlevel
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmsetapplicationcontext_ DMSETAPPLICATIONCONTEXT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmsetapplicationcontext_ dmsetapplicationcontext
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetapplicationcontext_ DMGETAPPLICATIONCONTEXT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetapplicationcontext_ dmgetapplicationcontext
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmhasvariablebounds_ DMHASVARIABLEBOUNDS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmhasvariablebounds_ dmhasvariablebounds
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmhascoloring_ DMHASCOLORING
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmhascoloring_ dmhascoloring
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmhascreaterestriction_ DMHASCREATERESTRICTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmhascreaterestriction_ dmhascreaterestriction
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmhascreateinjection_ DMHASCREATEINJECTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmhascreateinjection_ dmhascreateinjection
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetlocalboundingbox_ DMGETLOCALBOUNDINGBOX
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetlocalboundingbox_ dmgetlocalboundingbox
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetboundingbox_ DMGETBOUNDINGBOX
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetboundingbox_ dmgetboundingbox
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetsection_ DMGETSECTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetsection_ dmgetsection
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetlocalsection_ DMGETLOCALSECTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetlocalsection_ dmgetlocalsection
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmsetsection_ DMSETSECTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmsetsection_ dmsetsection
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmsetlocalsection_ DMSETLOCALSECTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmsetlocalsection_ dmsetlocalsection
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetdefaultconstraints_ DMGETDEFAULTCONSTRAINTS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetdefaultconstraints_ dmgetdefaultconstraints
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmsetdefaultconstraints_ DMSETDEFAULTCONSTRAINTS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmsetdefaultconstraints_ dmsetdefaultconstraints
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetglobalsection_ DMGETGLOBALSECTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetglobalsection_ dmgetglobalsection
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmsetglobalsection_ DMSETGLOBALSECTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmsetglobalsection_ dmsetglobalsection
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetsectionsf_ DMGETSECTIONSF
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetsectionsf_ dmgetsectionsf
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmsetsectionsf_ DMSETSECTIONSF
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmsetsectionsf_ dmsetsectionsf
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetpointsf_ DMGETPOINTSF
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetpointsf_ dmgetpointsf
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmsetpointsf_ DMSETPOINTSF
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmsetpointsf_ dmsetpointsf
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmclearfields_ DMCLEARFIELDS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmclearfields_ dmclearfields
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetnumfields_ DMGETNUMFIELDS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetnumfields_ dmgetnumfields
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmsetnumfields_ DMSETNUMFIELDS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmsetnumfields_ dmsetnumfields
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetfield_ DMGETFIELD
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetfield_ dmgetfield
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmsetfield_ DMSETFIELD
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmsetfield_ dmsetfield
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmaddfield_ DMADDFIELD
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmaddfield_ dmaddfield
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmsetfieldavoidtensor_ DMSETFIELDAVOIDTENSOR
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmsetfieldavoidtensor_ dmsetfieldavoidtensor
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetfieldavoidtensor_ DMGETFIELDAVOIDTENSOR
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetfieldavoidtensor_ dmgetfieldavoidtensor
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmcopyfields_ DMCOPYFIELDS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmcopyfields_ dmcopyfields
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetadjacency_ DMGETADJACENCY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetadjacency_ dmgetadjacency
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmsetadjacency_ DMSETADJACENCY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmsetadjacency_ dmsetadjacency
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetbasicadjacency_ DMGETBASICADJACENCY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetbasicadjacency_ dmgetbasicadjacency
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmsetbasicadjacency_ DMSETBASICADJACENCY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmsetbasicadjacency_ dmsetbasicadjacency
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetnumds_ DMGETNUMDS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetnumds_ dmgetnumds
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmcleards_ DMCLEARDS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmcleards_ dmcleards
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetds_ DMGETDS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetds_ dmgetds
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetcellds_ DMGETCELLDS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetcellds_ dmgetcellds
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetregionds_ DMGETREGIONDS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetregionds_ dmgetregionds
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmsetregionds_ DMSETREGIONDS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmsetregionds_ dmsetregionds
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetregionnumds_ DMGETREGIONNUMDS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetregionnumds_ dmgetregionnumds
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmsetregionnumds_ DMSETREGIONNUMDS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmsetregionnumds_ dmsetregionnumds
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmfindregionnum_ DMFINDREGIONNUM
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmfindregionnum_ dmfindregionnum
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmcreateds_ DMCREATEDS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmcreateds_ dmcreateds
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmcomputeexactsolution_ DMCOMPUTEEXACTSOLUTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmcomputeexactsolution_ dmcomputeexactsolution
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmcopyds_ DMCOPYDS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmcopyds_ dmcopyds
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmcopydisc_ DMCOPYDISC
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmcopydisc_ dmcopydisc
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetdimension_ DMGETDIMENSION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetdimension_ dmgetdimension
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmsetdimension_ DMSETDIMENSION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmsetdimension_ dmsetdimension
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetdimpoints_ DMGETDIMPOINTS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetdimpoints_ dmgetdimpoints
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmsetcoordinates_ DMSETCOORDINATES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmsetcoordinates_ dmsetcoordinates
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmsetcoordinateslocal_ DMSETCOORDINATESLOCAL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmsetcoordinateslocal_ dmsetcoordinateslocal
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetcoordinates_ DMGETCOORDINATES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetcoordinates_ dmgetcoordinates
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetcoordinateslocalsetup_ DMGETCOORDINATESLOCALSETUP
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetcoordinateslocalsetup_ dmgetcoordinateslocalsetup
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetcoordinateslocal_ DMGETCOORDINATESLOCAL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetcoordinateslocal_ dmgetcoordinateslocal
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetcoordinateslocalnoncollective_ DMGETCOORDINATESLOCALNONCOLLECTIVE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetcoordinateslocalnoncollective_ dmgetcoordinateslocalnoncollective
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetcoordinateslocaltuple_ DMGETCOORDINATESLOCALTUPLE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetcoordinateslocaltuple_ dmgetcoordinateslocaltuple
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetcoordinatedm_ DMGETCOORDINATEDM
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetcoordinatedm_ dmgetcoordinatedm
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmsetcoordinatedm_ DMSETCOORDINATEDM
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmsetcoordinatedm_ dmsetcoordinatedm
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetcoordinatedim_ DMGETCOORDINATEDIM
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetcoordinatedim_ dmgetcoordinatedim
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmsetcoordinatedim_ DMSETCOORDINATEDIM
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmsetcoordinatedim_ dmsetcoordinatedim
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetcoordinatesection_ DMGETCOORDINATESECTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetcoordinatesection_ dmgetcoordinatesection
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmsetcoordinatesection_ DMSETCOORDINATESECTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmsetcoordinatesection_ dmsetcoordinatesection
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmprojectcoordinates_ DMPROJECTCOORDINATES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmprojectcoordinates_ dmprojectcoordinates
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmlocalizecoordinate_ DMLOCALIZECOORDINATE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmlocalizecoordinate_ dmlocalizecoordinate
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetcoordinateslocalizedlocal_ DMGETCOORDINATESLOCALIZEDLOCAL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetcoordinateslocalizedlocal_ dmgetcoordinateslocalizedlocal
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetcoordinateslocalized_ DMGETCOORDINATESLOCALIZED
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetcoordinateslocalized_ dmgetcoordinateslocalized
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmlocalizecoordinates_ DMLOCALIZECOORDINATES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmlocalizecoordinates_ dmlocalizecoordinates
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmlocatepoints_ DMLOCATEPOINTS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmlocatepoints_ dmlocatepoints
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetoutputdm_ DMGETOUTPUTDM
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetoutputdm_ dmgetoutputdm
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetoutputsequencenumber_ DMGETOUTPUTSEQUENCENUMBER
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetoutputsequencenumber_ dmgetoutputsequencenumber
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmsetoutputsequencenumber_ DMSETOUTPUTSEQUENCENUMBER
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmsetoutputsequencenumber_ dmsetoutputsequencenumber
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetusenatural_ DMGETUSENATURAL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetusenatural_ dmgetusenatural
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmsetusenatural_ DMSETUSENATURAL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmsetusenatural_ dmsetusenatural
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetnumlabels_ DMGETNUMLABELS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetnumlabels_ dmgetnumlabels
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmremovelabelbyself_ DMREMOVELABELBYSELF
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmremovelabelbyself_ dmremovelabelbyself
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmcopylabels_ DMCOPYLABELS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmcopylabels_ dmcopylabels
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetcoarsedm_ DMGETCOARSEDM
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetcoarsedm_ dmgetcoarsedm
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmsetcoarsedm_ DMSETCOARSEDM
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmsetcoarsedm_ dmsetcoarsedm
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetfinedm_ DMGETFINEDM
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetfinedm_ dmgetfinedm
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmsetfinedm_ DMSETFINEDM
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmsetfinedm_ dmsetfinedm
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetnumboundary_ DMGETNUMBOUNDARY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetnumboundary_ dmgetnumboundary
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define matfdcoloringusedm_ MATFDCOLORINGUSEDM
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define matfdcoloringusedm_ matfdcoloringusedm
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmgetcompatibility_ DMGETCOMPATIBILITY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmgetcompatibility_ dmgetcompatibility
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmmonitorcancel_ DMMONITORCANCEL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmmonitorcancel_ dmmonitorcancel
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmmonitor_ DMMONITOR
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmmonitor_ dmmonitor
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmcomputeerror_ DMCOMPUTEERROR
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmcomputeerror_ dmcomputeerror
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
PETSC_EXTERN void  dmcreate_(MPI_Fint * comm,DM *dm, int *__ierr)
{
*__ierr = DMCreate(
	MPI_Comm_f2c(*(comm)),dm);
}
PETSC_EXTERN void  dmclone_(DM dm,DM *newdm, int *__ierr)
{
*__ierr = DMClone(
	(DM)PetscToPointer((dm) ),newdm);
}
PETSC_EXTERN void  vecgetdm_(Vec v,DM *dm, int *__ierr)
{
*__ierr = VecGetDM(
	(Vec)PetscToPointer((v) ),dm);
}
PETSC_EXTERN void  vecsetdm_(Vec v,DM dm, int *__ierr)
{
*__ierr = VecSetDM(
	(Vec)PetscToPointer((v) ),
	(DM)PetscToPointer((dm) ));
}
PETSC_EXTERN void  matgetdm_(Mat A,DM *dm, int *__ierr)
{
*__ierr = MatGetDM(
	(Mat)PetscToPointer((A) ),dm);
}
PETSC_EXTERN void  matsetdm_(Mat A,DM dm, int *__ierr)
{
*__ierr = MatSetDM(
	(Mat)PetscToPointer((A) ),
	(DM)PetscToPointer((dm) ));
}
PETSC_EXTERN void  dmsetup_(DM dm, int *__ierr)
{
*__ierr = DMSetUp(
	(DM)PetscToPointer((dm) ));
}
PETSC_EXTERN void  dmsetfromoptions_(DM dm, int *__ierr)
{
*__ierr = DMSetFromOptions(
	(DM)PetscToPointer((dm) ));
}
PETSC_EXTERN void  dmcreateglobalvector_(DM dm,Vec *vec, int *__ierr)
{
*__ierr = DMCreateGlobalVector(
	(DM)PetscToPointer((dm) ),vec);
}
PETSC_EXTERN void  dmcreatelocalvector_(DM dm,Vec *vec, int *__ierr)
{
*__ierr = DMCreateLocalVector(
	(DM)PetscToPointer((dm) ),vec);
}
PETSC_EXTERN void  dmgetlocaltoglobalmapping_(DM dm,ISLocalToGlobalMapping *ltog, int *__ierr)
{
*__ierr = DMGetLocalToGlobalMapping(
	(DM)PetscToPointer((dm) ),ltog);
}
PETSC_EXTERN void  dmgetblocksize_(DM dm,PetscInt *bs, int *__ierr)
{
*__ierr = DMGetBlockSize(
	(DM)PetscToPointer((dm) ),bs);
}
PETSC_EXTERN void  dmcreateinterpolationscale_(DM dac,DM daf,Mat mat,Vec *scale, int *__ierr)
{
*__ierr = DMCreateInterpolationScale(
	(DM)PetscToPointer((dac) ),
	(DM)PetscToPointer((daf) ),
	(Mat)PetscToPointer((mat) ),scale);
}
PETSC_EXTERN void  dmcreaterestriction_(DM dmc,DM dmf,Mat *mat, int *__ierr)
{
*__ierr = DMCreateRestriction(
	(DM)PetscToPointer((dmc) ),
	(DM)PetscToPointer((dmf) ),mat);
}
PETSC_EXTERN void  dmcreateinjection_(DM dac,DM daf,Mat *mat, int *__ierr)
{
*__ierr = DMCreateInjection(
	(DM)PetscToPointer((dac) ),
	(DM)PetscToPointer((daf) ),mat);
}
PETSC_EXTERN void  dmcreatemassmatrix_(DM dac,DM daf,Mat *mat, int *__ierr)
{
*__ierr = DMCreateMassMatrix(
	(DM)PetscToPointer((dac) ),
	(DM)PetscToPointer((daf) ),mat);
}
PETSC_EXTERN void  dmcreatecoloring_(DM dm,ISColoringType *ctype,ISColoring *coloring, int *__ierr)
{
*__ierr = DMCreateColoring(
	(DM)PetscToPointer((dm) ),*ctype,coloring);
}
PETSC_EXTERN void  dmcreatematrix_(DM dm,Mat *mat, int *__ierr)
{
*__ierr = DMCreateMatrix(
	(DM)PetscToPointer((dm) ),mat);
}
PETSC_EXTERN void  dmsetmatrixpreallocateonly_(DM dm,PetscBool *only, int *__ierr)
{
*__ierr = DMSetMatrixPreallocateOnly(
	(DM)PetscToPointer((dm) ),*only);
}
PETSC_EXTERN void  dmsetmatrixstructureonly_(DM dm,PetscBool *only, int *__ierr)
{
*__ierr = DMSetMatrixStructureOnly(
	(DM)PetscToPointer((dm) ),*only);
}
PETSC_EXTERN void  dmcreatesubdm_(DM dm,PetscInt *numFields, PetscInt fields[],IS *is,DM *subdm, int *__ierr)
{
*__ierr = DMCreateSubDM(
	(DM)PetscToPointer((dm) ),*numFields,fields,is,subdm);
}
PETSC_EXTERN void  dmrefine_(DM dm,MPI_Fint * comm,DM *dmf, int *__ierr)
{
*__ierr = DMRefine(
	(DM)PetscToPointer((dm) ),
	MPI_Comm_f2c(*(comm)),dmf);
}
PETSC_EXTERN void  dminterpolate_(DM coarse,Mat interp,DM fine, int *__ierr)
{
*__ierr = DMInterpolate(
	(DM)PetscToPointer((coarse) ),
	(Mat)PetscToPointer((interp) ),
	(DM)PetscToPointer((fine) ));
}
PETSC_EXTERN void  dminterpolatesolution_(DM coarse,DM fine,Mat interp,Vec coarseSol,Vec fineSol, int *__ierr)
{
*__ierr = DMInterpolateSolution(
	(DM)PetscToPointer((coarse) ),
	(DM)PetscToPointer((fine) ),
	(Mat)PetscToPointer((interp) ),
	(Vec)PetscToPointer((coarseSol) ),
	(Vec)PetscToPointer((fineSol) ));
}
PETSC_EXTERN void  dmgetrefinelevel_(DM dm,PetscInt *level, int *__ierr)
{
*__ierr = DMGetRefineLevel(
	(DM)PetscToPointer((dm) ),level);
}
PETSC_EXTERN void  dmsetrefinelevel_(DM dm,PetscInt *level, int *__ierr)
{
*__ierr = DMSetRefineLevel(
	(DM)PetscToPointer((dm) ),*level);
}
PETSC_EXTERN void  dmhasbasistransform_(DM dm,PetscBool *flg, int *__ierr)
{
*__ierr = DMHasBasisTransform(
	(DM)PetscToPointer((dm) ),flg);
}
PETSC_EXTERN void  dmglobaltolocal_(DM dm,Vec g,InsertMode *mode,Vec l, int *__ierr)
{
*__ierr = DMGlobalToLocal(
	(DM)PetscToPointer((dm) ),
	(Vec)PetscToPointer((g) ),*mode,
	(Vec)PetscToPointer((l) ));
}
PETSC_EXTERN void  dmglobaltolocalbegin_(DM dm,Vec g,InsertMode *mode,Vec l, int *__ierr)
{
*__ierr = DMGlobalToLocalBegin(
	(DM)PetscToPointer((dm) ),
	(Vec)PetscToPointer((g) ),*mode,
	(Vec)PetscToPointer((l) ));
}
PETSC_EXTERN void  dmglobaltolocalend_(DM dm,Vec g,InsertMode *mode,Vec l, int *__ierr)
{
*__ierr = DMGlobalToLocalEnd(
	(DM)PetscToPointer((dm) ),
	(Vec)PetscToPointer((g) ),*mode,
	(Vec)PetscToPointer((l) ));
}
PETSC_EXTERN void  dmlocaltoglobal_(DM dm,Vec l,InsertMode *mode,Vec g, int *__ierr)
{
*__ierr = DMLocalToGlobal(
	(DM)PetscToPointer((dm) ),
	(Vec)PetscToPointer((l) ),*mode,
	(Vec)PetscToPointer((g) ));
}
PETSC_EXTERN void  dmlocaltoglobalbegin_(DM dm,Vec l,InsertMode *mode,Vec g, int *__ierr)
{
*__ierr = DMLocalToGlobalBegin(
	(DM)PetscToPointer((dm) ),
	(Vec)PetscToPointer((l) ),*mode,
	(Vec)PetscToPointer((g) ));
}
PETSC_EXTERN void  dmlocaltoglobalend_(DM dm,Vec l,InsertMode *mode,Vec g, int *__ierr)
{
*__ierr = DMLocalToGlobalEnd(
	(DM)PetscToPointer((dm) ),
	(Vec)PetscToPointer((l) ),*mode,
	(Vec)PetscToPointer((g) ));
}
PETSC_EXTERN void  dmlocaltolocalbegin_(DM dm,Vec g,InsertMode *mode,Vec l, int *__ierr)
{
*__ierr = DMLocalToLocalBegin(
	(DM)PetscToPointer((dm) ),
	(Vec)PetscToPointer((g) ),*mode,
	(Vec)PetscToPointer((l) ));
}
PETSC_EXTERN void  dmlocaltolocalend_(DM dm,Vec g,InsertMode *mode,Vec l, int *__ierr)
{
*__ierr = DMLocalToLocalEnd(
	(DM)PetscToPointer((dm) ),
	(Vec)PetscToPointer((g) ),*mode,
	(Vec)PetscToPointer((l) ));
}
PETSC_EXTERN void  dmcoarsen_(DM dm,MPI_Fint * comm,DM *dmc, int *__ierr)
{
*__ierr = DMCoarsen(
	(DM)PetscToPointer((dm) ),
	MPI_Comm_f2c(*(comm)),dmc);
}
PETSC_EXTERN void  dmrestrict_(DM fine,Mat restrct,Vec rscale,Mat inject,DM coarse, int *__ierr)
{
*__ierr = DMRestrict(
	(DM)PetscToPointer((fine) ),
	(Mat)PetscToPointer((restrct) ),
	(Vec)PetscToPointer((rscale) ),
	(Mat)PetscToPointer((inject) ),
	(DM)PetscToPointer((coarse) ));
}
PETSC_EXTERN void  dmsubdomainrestrict_(DM global,VecScatter oscatter,VecScatter gscatter,DM subdm, int *__ierr)
{
*__ierr = DMSubDomainRestrict(
	(DM)PetscToPointer((global) ),
	(VecScatter)PetscToPointer((oscatter) ),
	(VecScatter)PetscToPointer((gscatter) ),
	(DM)PetscToPointer((subdm) ));
}
PETSC_EXTERN void  dmgetcoarsenlevel_(DM dm,PetscInt *level, int *__ierr)
{
*__ierr = DMGetCoarsenLevel(
	(DM)PetscToPointer((dm) ),level);
}
PETSC_EXTERN void  dmsetcoarsenlevel_(DM dm,PetscInt *level, int *__ierr)
{
*__ierr = DMSetCoarsenLevel(
	(DM)PetscToPointer((dm) ),*level);
}
PETSC_EXTERN void  dmsetapplicationcontext_(DM dm,void*ctx, int *__ierr)
{
*__ierr = DMSetApplicationContext(
	(DM)PetscToPointer((dm) ),ctx);
}
PETSC_EXTERN void  dmgetapplicationcontext_(DM dm,void*ctx, int *__ierr)
{
*__ierr = DMGetApplicationContext(
	(DM)PetscToPointer((dm) ),ctx);
}
PETSC_EXTERN void  dmhasvariablebounds_(DM dm,PetscBool *flg, int *__ierr)
{
*__ierr = DMHasVariableBounds(
	(DM)PetscToPointer((dm) ),flg);
}
PETSC_EXTERN void  dmhascoloring_(DM dm,PetscBool *flg, int *__ierr)
{
*__ierr = DMHasColoring(
	(DM)PetscToPointer((dm) ),flg);
}
PETSC_EXTERN void  dmhascreaterestriction_(DM dm,PetscBool *flg, int *__ierr)
{
*__ierr = DMHasCreateRestriction(
	(DM)PetscToPointer((dm) ),flg);
}
PETSC_EXTERN void  dmhascreateinjection_(DM dm,PetscBool *flg, int *__ierr)
{
*__ierr = DMHasCreateInjection(
	(DM)PetscToPointer((dm) ),flg);
}
PETSC_EXTERN void  dmgetlocalboundingbox_(DM dm,PetscReal lmin[],PetscReal lmax[], int *__ierr)
{
*__ierr = DMGetLocalBoundingBox(
	(DM)PetscToPointer((dm) ),lmin,lmax);
}
PETSC_EXTERN void  dmgetboundingbox_(DM dm,PetscReal gmin[],PetscReal gmax[], int *__ierr)
{
*__ierr = DMGetBoundingBox(
	(DM)PetscToPointer((dm) ),gmin,gmax);
}
PETSC_EXTERN void  dmgetsection_(DM dm,PetscSection *section, int *__ierr)
{
*__ierr = DMGetSection(
	(DM)PetscToPointer((dm) ),section);
}
PETSC_EXTERN void  dmgetlocalsection_(DM dm,PetscSection *section, int *__ierr)
{
*__ierr = DMGetLocalSection(
	(DM)PetscToPointer((dm) ),section);
}
PETSC_EXTERN void  dmsetsection_(DM dm,PetscSection section, int *__ierr)
{
*__ierr = DMSetSection(
	(DM)PetscToPointer((dm) ),
	(PetscSection)PetscToPointer((section) ));
}
PETSC_EXTERN void  dmsetlocalsection_(DM dm,PetscSection section, int *__ierr)
{
*__ierr = DMSetLocalSection(
	(DM)PetscToPointer((dm) ),
	(PetscSection)PetscToPointer((section) ));
}
PETSC_EXTERN void  dmgetdefaultconstraints_(DM dm,PetscSection *section,Mat *mat, int *__ierr)
{
*__ierr = DMGetDefaultConstraints(
	(DM)PetscToPointer((dm) ),section,mat);
}
PETSC_EXTERN void  dmsetdefaultconstraints_(DM dm,PetscSection section,Mat mat, int *__ierr)
{
*__ierr = DMSetDefaultConstraints(
	(DM)PetscToPointer((dm) ),
	(PetscSection)PetscToPointer((section) ),
	(Mat)PetscToPointer((mat) ));
}
PETSC_EXTERN void  dmgetglobalsection_(DM dm,PetscSection *section, int *__ierr)
{
*__ierr = DMGetGlobalSection(
	(DM)PetscToPointer((dm) ),section);
}
PETSC_EXTERN void  dmsetglobalsection_(DM dm,PetscSection section, int *__ierr)
{
*__ierr = DMSetGlobalSection(
	(DM)PetscToPointer((dm) ),
	(PetscSection)PetscToPointer((section) ));
}
PETSC_EXTERN void  dmgetsectionsf_(DM dm,PetscSF *sf, int *__ierr)
{
*__ierr = DMGetSectionSF(
	(DM)PetscToPointer((dm) ),sf);
}
PETSC_EXTERN void  dmsetsectionsf_(DM dm,PetscSF sf, int *__ierr)
{
*__ierr = DMSetSectionSF(
	(DM)PetscToPointer((dm) ),
	(PetscSF)PetscToPointer((sf) ));
}
PETSC_EXTERN void  dmgetpointsf_(DM dm,PetscSF *sf, int *__ierr)
{
*__ierr = DMGetPointSF(
	(DM)PetscToPointer((dm) ),sf);
}
PETSC_EXTERN void  dmsetpointsf_(DM dm,PetscSF sf, int *__ierr)
{
*__ierr = DMSetPointSF(
	(DM)PetscToPointer((dm) ),
	(PetscSF)PetscToPointer((sf) ));
}
PETSC_EXTERN void  dmclearfields_(DM dm, int *__ierr)
{
*__ierr = DMClearFields(
	(DM)PetscToPointer((dm) ));
}
PETSC_EXTERN void  dmgetnumfields_(DM dm,PetscInt *numFields, int *__ierr)
{
*__ierr = DMGetNumFields(
	(DM)PetscToPointer((dm) ),numFields);
}
PETSC_EXTERN void  dmsetnumfields_(DM dm,PetscInt *numFields, int *__ierr)
{
*__ierr = DMSetNumFields(
	(DM)PetscToPointer((dm) ),*numFields);
}
PETSC_EXTERN void  dmgetfield_(DM dm,PetscInt *f,DMLabel *label,PetscObject *field, int *__ierr)
{
*__ierr = DMGetField(
	(DM)PetscToPointer((dm) ),*f,label,field);
}
PETSC_EXTERN void  dmsetfield_(DM dm,PetscInt *f,DMLabel label,PetscObject field, int *__ierr)
{
*__ierr = DMSetField(
	(DM)PetscToPointer((dm) ),*f,
	(DMLabel)PetscToPointer((label) ),
	(PetscObject)PetscToPointer((field) ));
}
PETSC_EXTERN void  dmaddfield_(DM dm,DMLabel label,PetscObject field, int *__ierr)
{
*__ierr = DMAddField(
	(DM)PetscToPointer((dm) ),
	(DMLabel)PetscToPointer((label) ),
	(PetscObject)PetscToPointer((field) ));
}
PETSC_EXTERN void  dmsetfieldavoidtensor_(DM dm,PetscInt *f,PetscBool *avoidTensor, int *__ierr)
{
*__ierr = DMSetFieldAvoidTensor(
	(DM)PetscToPointer((dm) ),*f,*avoidTensor);
}
PETSC_EXTERN void  dmgetfieldavoidtensor_(DM dm,PetscInt *f,PetscBool *avoidTensor, int *__ierr)
{
*__ierr = DMGetFieldAvoidTensor(
	(DM)PetscToPointer((dm) ),*f,avoidTensor);
}
PETSC_EXTERN void  dmcopyfields_(DM dm,DM newdm, int *__ierr)
{
*__ierr = DMCopyFields(
	(DM)PetscToPointer((dm) ),
	(DM)PetscToPointer((newdm) ));
}
PETSC_EXTERN void  dmgetadjacency_(DM dm,PetscInt *f,PetscBool *useCone,PetscBool *useClosure, int *__ierr)
{
*__ierr = DMGetAdjacency(
	(DM)PetscToPointer((dm) ),*f,useCone,useClosure);
}
PETSC_EXTERN void  dmsetadjacency_(DM dm,PetscInt *f,PetscBool *useCone,PetscBool *useClosure, int *__ierr)
{
*__ierr = DMSetAdjacency(
	(DM)PetscToPointer((dm) ),*f,*useCone,*useClosure);
}
PETSC_EXTERN void  dmgetbasicadjacency_(DM dm,PetscBool *useCone,PetscBool *useClosure, int *__ierr)
{
*__ierr = DMGetBasicAdjacency(
	(DM)PetscToPointer((dm) ),useCone,useClosure);
}
PETSC_EXTERN void  dmsetbasicadjacency_(DM dm,PetscBool *useCone,PetscBool *useClosure, int *__ierr)
{
*__ierr = DMSetBasicAdjacency(
	(DM)PetscToPointer((dm) ),*useCone,*useClosure);
}
PETSC_EXTERN void  dmgetnumds_(DM dm,PetscInt *Nds, int *__ierr)
{
*__ierr = DMGetNumDS(
	(DM)PetscToPointer((dm) ),Nds);
}
PETSC_EXTERN void  dmcleards_(DM dm, int *__ierr)
{
*__ierr = DMClearDS(
	(DM)PetscToPointer((dm) ));
}
PETSC_EXTERN void  dmgetds_(DM dm,PetscDS *prob, int *__ierr)
{
*__ierr = DMGetDS(
	(DM)PetscToPointer((dm) ),prob);
}
PETSC_EXTERN void  dmgetcellds_(DM dm,PetscInt *point,PetscDS *prob, int *__ierr)
{
*__ierr = DMGetCellDS(
	(DM)PetscToPointer((dm) ),*point,prob);
}
PETSC_EXTERN void  dmgetregionds_(DM dm,DMLabel label,IS *fields,PetscDS *ds, int *__ierr)
{
*__ierr = DMGetRegionDS(
	(DM)PetscToPointer((dm) ),
	(DMLabel)PetscToPointer((label) ),fields,ds);
}
PETSC_EXTERN void  dmsetregionds_(DM dm,DMLabel label,IS fields,PetscDS ds, int *__ierr)
{
*__ierr = DMSetRegionDS(
	(DM)PetscToPointer((dm) ),
	(DMLabel)PetscToPointer((label) ),
	(IS)PetscToPointer((fields) ),
	(PetscDS)PetscToPointer((ds) ));
}
PETSC_EXTERN void  dmgetregionnumds_(DM dm,PetscInt *num,DMLabel *label,IS *fields,PetscDS *ds, int *__ierr)
{
*__ierr = DMGetRegionNumDS(
	(DM)PetscToPointer((dm) ),*num,label,fields,ds);
}
PETSC_EXTERN void  dmsetregionnumds_(DM dm,PetscInt *num,DMLabel label,IS fields,PetscDS ds, int *__ierr)
{
*__ierr = DMSetRegionNumDS(
	(DM)PetscToPointer((dm) ),*num,
	(DMLabel)PetscToPointer((label) ),
	(IS)PetscToPointer((fields) ),
	(PetscDS)PetscToPointer((ds) ));
}
PETSC_EXTERN void  dmfindregionnum_(DM dm,PetscDS ds,PetscInt *num, int *__ierr)
{
*__ierr = DMFindRegionNum(
	(DM)PetscToPointer((dm) ),
	(PetscDS)PetscToPointer((ds) ),num);
}
PETSC_EXTERN void  dmcreateds_(DM dm, int *__ierr)
{
*__ierr = DMCreateDS(
	(DM)PetscToPointer((dm) ));
}
PETSC_EXTERN void  dmcomputeexactsolution_(DM dm,PetscReal *time,Vec u,Vec u_t, int *__ierr)
{
*__ierr = DMComputeExactSolution(
	(DM)PetscToPointer((dm) ),*time,
	(Vec)PetscToPointer((u) ),
	(Vec)PetscToPointer((u_t) ));
}
PETSC_EXTERN void  dmcopyds_(DM dm,DM newdm, int *__ierr)
{
*__ierr = DMCopyDS(
	(DM)PetscToPointer((dm) ),
	(DM)PetscToPointer((newdm) ));
}
PETSC_EXTERN void  dmcopydisc_(DM dm,DM newdm, int *__ierr)
{
*__ierr = DMCopyDisc(
	(DM)PetscToPointer((dm) ),
	(DM)PetscToPointer((newdm) ));
}
PETSC_EXTERN void  dmgetdimension_(DM dm,PetscInt *dim, int *__ierr)
{
*__ierr = DMGetDimension(
	(DM)PetscToPointer((dm) ),dim);
}
PETSC_EXTERN void  dmsetdimension_(DM dm,PetscInt *dim, int *__ierr)
{
*__ierr = DMSetDimension(
	(DM)PetscToPointer((dm) ),*dim);
}
PETSC_EXTERN void  dmgetdimpoints_(DM dm,PetscInt *dim,PetscInt *pStart,PetscInt *pEnd, int *__ierr)
{
*__ierr = DMGetDimPoints(
	(DM)PetscToPointer((dm) ),*dim,pStart,pEnd);
}
PETSC_EXTERN void  dmsetcoordinates_(DM dm,Vec c, int *__ierr)
{
*__ierr = DMSetCoordinates(
	(DM)PetscToPointer((dm) ),
	(Vec)PetscToPointer((c) ));
}
PETSC_EXTERN void  dmsetcoordinateslocal_(DM dm,Vec c, int *__ierr)
{
*__ierr = DMSetCoordinatesLocal(
	(DM)PetscToPointer((dm) ),
	(Vec)PetscToPointer((c) ));
}
PETSC_EXTERN void  dmgetcoordinates_(DM dm,Vec *c, int *__ierr)
{
*__ierr = DMGetCoordinates(
	(DM)PetscToPointer((dm) ),c);
}
PETSC_EXTERN void  dmgetcoordinateslocalsetup_(DM dm, int *__ierr)
{
*__ierr = DMGetCoordinatesLocalSetUp(
	(DM)PetscToPointer((dm) ));
}
PETSC_EXTERN void  dmgetcoordinateslocal_(DM dm,Vec *c, int *__ierr)
{
*__ierr = DMGetCoordinatesLocal(
	(DM)PetscToPointer((dm) ),c);
}
PETSC_EXTERN void  dmgetcoordinateslocalnoncollective_(DM dm,Vec *c, int *__ierr)
{
*__ierr = DMGetCoordinatesLocalNoncollective(
	(DM)PetscToPointer((dm) ),c);
}
PETSC_EXTERN void  dmgetcoordinateslocaltuple_(DM dm,IS p,PetscSection *pCoordSection,Vec *pCoord, int *__ierr)
{
*__ierr = DMGetCoordinatesLocalTuple(
	(DM)PetscToPointer((dm) ),
	(IS)PetscToPointer((p) ),pCoordSection,pCoord);
}
PETSC_EXTERN void  dmgetcoordinatedm_(DM dm,DM *cdm, int *__ierr)
{
*__ierr = DMGetCoordinateDM(
	(DM)PetscToPointer((dm) ),cdm);
}
PETSC_EXTERN void  dmsetcoordinatedm_(DM dm,DM cdm, int *__ierr)
{
*__ierr = DMSetCoordinateDM(
	(DM)PetscToPointer((dm) ),
	(DM)PetscToPointer((cdm) ));
}
PETSC_EXTERN void  dmgetcoordinatedim_(DM dm,PetscInt *dim, int *__ierr)
{
*__ierr = DMGetCoordinateDim(
	(DM)PetscToPointer((dm) ),dim);
}
PETSC_EXTERN void  dmsetcoordinatedim_(DM dm,PetscInt *dim, int *__ierr)
{
*__ierr = DMSetCoordinateDim(
	(DM)PetscToPointer((dm) ),*dim);
}
PETSC_EXTERN void  dmgetcoordinatesection_(DM dm,PetscSection *section, int *__ierr)
{
*__ierr = DMGetCoordinateSection(
	(DM)PetscToPointer((dm) ),section);
}
PETSC_EXTERN void  dmsetcoordinatesection_(DM dm,PetscInt *dim,PetscSection section, int *__ierr)
{
*__ierr = DMSetCoordinateSection(
	(DM)PetscToPointer((dm) ),*dim,
	(PetscSection)PetscToPointer((section) ));
}
PETSC_EXTERN void  dmprojectcoordinates_(DM dm,PetscFE disc, int *__ierr)
{
*__ierr = DMProjectCoordinates(
	(DM)PetscToPointer((dm) ),
	(PetscFE)PetscToPointer((disc) ));
}
PETSC_EXTERN void  dmlocalizecoordinate_(DM dm, PetscScalar in[],PetscBool *endpoint,PetscScalar out[], int *__ierr)
{
*__ierr = DMLocalizeCoordinate(
	(DM)PetscToPointer((dm) ),in,*endpoint,out);
}
PETSC_EXTERN void  dmgetcoordinateslocalizedlocal_(DM dm,PetscBool *areLocalized, int *__ierr)
{
*__ierr = DMGetCoordinatesLocalizedLocal(
	(DM)PetscToPointer((dm) ),areLocalized);
}
PETSC_EXTERN void  dmgetcoordinateslocalized_(DM dm,PetscBool *areLocalized, int *__ierr)
{
*__ierr = DMGetCoordinatesLocalized(
	(DM)PetscToPointer((dm) ),areLocalized);
}
PETSC_EXTERN void  dmlocalizecoordinates_(DM dm, int *__ierr)
{
*__ierr = DMLocalizeCoordinates(
	(DM)PetscToPointer((dm) ));
}
PETSC_EXTERN void  dmlocatepoints_(DM dm,Vec v,DMPointLocationType *ltype,PetscSF *cellSF, int *__ierr)
{
*__ierr = DMLocatePoints(
	(DM)PetscToPointer((dm) ),
	(Vec)PetscToPointer((v) ),*ltype,cellSF);
}
PETSC_EXTERN void  dmgetoutputdm_(DM dm,DM *odm, int *__ierr)
{
*__ierr = DMGetOutputDM(
	(DM)PetscToPointer((dm) ),odm);
}
PETSC_EXTERN void  dmgetoutputsequencenumber_(DM dm,PetscInt *num,PetscReal *val, int *__ierr)
{
*__ierr = DMGetOutputSequenceNumber(
	(DM)PetscToPointer((dm) ),num,val);
}
PETSC_EXTERN void  dmsetoutputsequencenumber_(DM dm,PetscInt *num,PetscReal *val, int *__ierr)
{
*__ierr = DMSetOutputSequenceNumber(
	(DM)PetscToPointer((dm) ),*num,*val);
}
PETSC_EXTERN void  dmgetusenatural_(DM dm,PetscBool *useNatural, int *__ierr)
{
*__ierr = DMGetUseNatural(
	(DM)PetscToPointer((dm) ),useNatural);
}
PETSC_EXTERN void  dmsetusenatural_(DM dm,PetscBool *useNatural, int *__ierr)
{
*__ierr = DMSetUseNatural(
	(DM)PetscToPointer((dm) ),*useNatural);
}
PETSC_EXTERN void  dmgetnumlabels_(DM dm,PetscInt *numLabels, int *__ierr)
{
*__ierr = DMGetNumLabels(
	(DM)PetscToPointer((dm) ),numLabels);
}
PETSC_EXTERN void  dmremovelabelbyself_(DM dm,DMLabel *label,PetscBool *failNotFound, int *__ierr)
{
*__ierr = DMRemoveLabelBySelf(
	(DM)PetscToPointer((dm) ),label,*failNotFound);
}
PETSC_EXTERN void  dmcopylabels_(DM dmA,DM dmB,PetscCopyMode *mode,PetscBool *all, int *__ierr)
{
*__ierr = DMCopyLabels(
	(DM)PetscToPointer((dmA) ),
	(DM)PetscToPointer((dmB) ),*mode,*all);
}
PETSC_EXTERN void  dmgetcoarsedm_(DM dm,DM *cdm, int *__ierr)
{
*__ierr = DMGetCoarseDM(
	(DM)PetscToPointer((dm) ),cdm);
}
PETSC_EXTERN void  dmsetcoarsedm_(DM dm,DM cdm, int *__ierr)
{
*__ierr = DMSetCoarseDM(
	(DM)PetscToPointer((dm) ),
	(DM)PetscToPointer((cdm) ));
}
PETSC_EXTERN void  dmgetfinedm_(DM dm,DM *fdm, int *__ierr)
{
*__ierr = DMGetFineDM(
	(DM)PetscToPointer((dm) ),fdm);
}
PETSC_EXTERN void  dmsetfinedm_(DM dm,DM fdm, int *__ierr)
{
*__ierr = DMSetFineDM(
	(DM)PetscToPointer((dm) ),
	(DM)PetscToPointer((fdm) ));
}
PETSC_EXTERN void  dmgetnumboundary_(DM dm,PetscInt *numBd, int *__ierr)
{
*__ierr = DMGetNumBoundary(
	(DM)PetscToPointer((dm) ),numBd);
}
PETSC_EXTERN void  matfdcoloringusedm_(Mat coloring,MatFDColoring fdcoloring, int *__ierr)
{
*__ierr = MatFDColoringUseDM(
	(Mat)PetscToPointer((coloring) ),
	(MatFDColoring)PetscToPointer((fdcoloring) ));
}

PETSC_EXTERN void  dmgetcompatibility_(DM dm1,DM dm2,PetscBool *compatible,PetscBool *set, int *__ierr)
{
*__ierr = DMGetCompatibility(
	(DM)PetscToPointer((dm1) ),
	(DM)PetscToPointer((dm2) ),compatible,set);
}
PETSC_EXTERN void  dmmonitorcancel_(DM dm, int *__ierr)
{
*__ierr = DMMonitorCancel(
	(DM)PetscToPointer((dm) ));
}
PETSC_EXTERN void  dmmonitor_(DM dm, int *__ierr)
{
*__ierr = DMMonitor(
	(DM)PetscToPointer((dm) ));
}
PETSC_EXTERN void  dmcomputeerror_(DM dm,Vec sol,PetscReal errors[],Vec *errorVec, int *__ierr)
{
*__ierr = DMComputeError(
	(DM)PetscToPointer((dm) ),
	(Vec)PetscToPointer((sol) ),errors,errorVec);
}
#if defined(__cplusplus)
}
#endif
