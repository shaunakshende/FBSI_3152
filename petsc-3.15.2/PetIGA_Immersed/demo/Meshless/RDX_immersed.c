#include "petiga.h"
#include "time.h"
#include <math.h>

typedef struct {
 PetscReal lamda,Lx,Ly,Lz,kappa,Cp,R,Cv,mu;
 PetscScalar F[4];
 IGA iga;
 IGA iga_energy;
 IGA iga_strain;
 PetscReal temp0;
 PetscReal p0;
 PetscReal rhoRDX;
 PetscReal currentTime;
 PetscReal timeStep;
 PetscReal finalTime;
 PetscInt  stepNumber;

 PetscReal Alpha_m,Alpha_f,Gamma,Beta;


 Vec A1,Ap,An,V1,Vp,Vn,D1,Dp,Dn,Aa,Va,Da,V0,A0,D0,dA;
 PetscInt  Nx,Ny,Nz,max_its;


 PetscReal TimeRestart;
 PetscInt  StepRestart;
 PetscInt  FreqRestarts;
 PetscInt  FreqResults;
 PetscReal H[3];
 PetscReal totalCurrentExplosiveVolume;
 PetscReal totalInitialExplosiveVolume;

} AppCtx;

//Information we need on the particles
typedef struct
{
	  PetscInt  		    nodeID;  //ID of the particle
	  PetscReal         nodalVolume; //nodalVolume associated with the particle. it is what we read from XVOL.dat and is used as the integration weigth when
								      //we perform nodal integration
    PetscReal         nodalDensity;
    PetscReal         nodalDensityInitial;
    PetscReal         nodalPressure;
    PetscReal         nodalVolumeInitial; //store the initial particle volume for update in explosive calculation
	  PetscReal         velocityGradient[4]; //velocity gradient on each particle. 2x2 for 2D. dvx/dx,dvx/dy,dvy/dx,dvy/dy for indexing 0-->4
    PetscReal         velocityGradientX[4];
    PetscReal 	  	  velocityGradientY[4];
    PetscReal         currentDeformationGradient[4];
    PetscReal         DeformationGradientOld[4];
    PetscReal         alphaDeformationGradient[4];
    PetscReal         determinantAlphaDeformationGradient;
} NODE;

//Information we need on the particles. Similar to the above. Everythin could be in one structure actually.
typedef struct
{
    int material; //0 for concrete, 1 for RDX
    PetscInt  ID; //ID of the particle. Same as nodeID above
    PetscReal initialCoord[2]; //Initial coordinates of the particle for 2D
    PetscReal currentCoord[2]; //Current coordinates of the particle for 2D
    PetscReal hvect[2];//dx, dy, and dz for each particle for 2D
    PetscReal totalPhysicalDisplacement[2]; //Displacement of the particle between two consecutive time steps in the two directions
    PetscReal totalPhysicalDisplacementOldIteration[2]; //Displacement of the particle in the previous iteration
    PetscReal totalPhysicalDisplacementOldStep[2]; //Displacement of the particle in the previous time step
    PetscReal totalPhysicalVelocity[2]; //Velocity of the particle at the current time position
    PetscReal totalPhysicalVelocityOldIteration[2]; //Velocity of the particle at the end of the previous iteration
    PetscReal totalPhysicalVelocityOldStep[2]; //Velocity of the particle at the end of the previous time step
    PetscReal totalPhysicalAcceleration[2]; //Acceleration of the particle at the current time position
    PetscReal totalPhysicalAccelerationOldIteration[2]; //Acceleration of the particle at the end of the previous iteration
    PetscReal totalPhysicalAccelerationOldStep[2]; //Acceleration of the particle at the end of the previous time step
    PetscReal AccelerationIncrement[2]; //Acceleration increment of the particle
    PetscReal totalStress[3]; //stress of particle at current time. For 2D we have 4 components but due to symmetry it collapses to 3 components.
    						  //The numbering is 0->sigmaxx, 1->sigmayy, 2->sigmaxy
    PetscReal totalStress0[3];//stress of particle at the end of the time step (after iterations have been completed and after we have rotated fully for objectivity)
    						  //Numbering same as above
    PetscReal totalStrain[3];//strain of the particle at current time. Same numbering as above
    PetscReal EpsilonPlastEq;//Equivalent plastic strain. Basically a non decreasing quantity that accumulates the plastic strain in every time step.
    PetscInt  Boundary;
    PetscReal totalStressderivativeX[3];
    PetscReal totalStressderivativeY[3];
    PetscReal totalStressderivativeX0[3];
    PetscReal totalStressderivativeY0[3];
    PetscReal pressure;
    PetscReal tempCoord[2];
    PetscReal U[2];
    PetscReal mapX[4];
    PetscReal mapU[4];
} POINTS;

typedef struct
{

	PetscReal initialTime; //Initial time of our computation. usually 0
	PetscReal finalTime;   //Final time of our computation. Problem dependent
	PetscReal currentTime; //Time currently (time in the time step we are in)
	PetscReal gamma; //parameter for time integration, not used for now
	PetscReal beta;  //parameter for time integration, not used for now
	PetscReal timeStep; //time step, problem and mesh dependent
        PetscInt  stepNumber; //number of the current step
	PetscReal sound_speed;
    PetscReal   youngModulus;
    PetscReal   poissonRatio;
    PetscReal   density;
    PetscReal   lambda; //Lame parameter
    PetscReal   mu;  //Lame parameter
    PetscReal   kappa;
    PetscReal   SigmaYinitial;  //The initial yield stress under which the material starts plasticizing. If we don't have hardening then this is constant
    							//throughout the computation
    PetscInt    FreqResults; //Frequency with wich we export results for post processing
    PetscInt    numPoints; //total number of particles
    PetscInt    numNodes;  //total number of particles
    PetscReal   Alpha_f;   //parameter for generalized aplha method

    POINTS      *puntos;  //structure POINT is a member of PARAMETERS
    NODE        *nodos;   //structure NODES is a member of PARAMETERS

} PARAMETERS;

#undef __FUNCT__
#define __FUNCT__ "ReadLastResults"
PetscErrorCode ReadLastResults(AppCtx *user, PARAMETERS *par,Vec U,Vec V,PetscInt StepRestart)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscInt i;
  MPI_Comm comm;
  PetscViewer viewer;

  // U
   char filename[256];
   sprintf(filename,"RestartU%d.dat",StepRestart);
   ierr = PetscObjectGetComm((PetscObject)U,&comm);CHKERRQ(ierr);
   ierr = PetscViewerBinaryOpen(comm,filename,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
   ierr = VecLoad(U,viewer);CHKERRQ(ierr);
   ierr = PetscViewerDestroy(&viewer);

   // V
   sprintf(filename,"RestartV%d.dat",StepRestart);
   ierr = PetscObjectGetComm((PetscObject)V,&comm);CHKERRQ(ierr);
   ierr = PetscViewerBinaryOpen(comm,filename,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
   ierr = VecLoad(V,viewer);CHKERRQ(ierr);
   ierr = PetscViewerDestroy(&viewer);

   char filenameTime[256];
   sprintf(filenameTime,"RestartTime%d.dat",StepRestart);
         if (  (fopen(filenameTime, "rt")) == NULL)  {
         PetscPrintf(PETSC_COMM_WORLD,"Time file NOT found!\n");
         } else {
         FILE *fL = fopen(filenameTime, "rt");
             ierr = fscanf(fL, "%le ",&par->currentTime);
         fclose(fL);
         }
         user->currentTime  = par->currentTime;

  	 // Velocity
  	  		char filenameVelocity[256];
  	  		sprintf(filenameVelocity,"RestartVelocity%d.dat",StepRestart);
  	  	  	    if (  (fopen(filenameVelocity, "rt")) == NULL)  {
  	  	  	    PetscPrintf(PETSC_COMM_WORLD,"Velocity file NOT found!\n");
  	  	  	    } else {
  	  	  	    FILE *fL = fopen(filenameVelocity, "rt");
  	  	          for (i=0;i<par->numNodes; i++){
  	  	        	  ierr = fscanf(fL, "%le %le",&par->puntos[i].totalPhysicalVelocity[0], &par->puntos[i].totalPhysicalVelocity[1]);
  	  	          }
  	  	  	    fclose(fL);
  	  	  	    }
  	  	  	    // Acceleration
  	  	  	    char filenameAcceleration[256];
	  	  	  		sprintf(filenameAcceleration,"RestartAcceleration%d.dat",StepRestart);

	  	  	  	  	    if (  (fopen(filenameAcceleration, "rt")) == NULL)  {
	  	  	  	  	    PetscPrintf(PETSC_COMM_WORLD,"Acceleration file NOT found!\n");
	  	  	  	  	    } else {

	  	  	  	  	    FILE *fL = fopen(filenameAcceleration, "rt");

	  	  	  	          for (i=0;i<par->numNodes; i++){
	  	  	  	        	  ierr = fscanf(fL, "%le %le",&par->puntos[i].totalPhysicalAcceleration[0], &par->puntos[i].totalPhysicalAcceleration[1]);
	  	  	  	          }
	  	  	  	  	    fclose(fL);
	  	  	  	  	    }

                      char filenameDensity[256];
                      sprintf(filenameDensity,"RestartDensity%d.dat",StepRestart);
                            if (  (fopen(filenameDensity, "rt")) == NULL)  {
                            PetscPrintf(PETSC_COMM_WORLD,"Density file NOT found!\n");
                            } else {
                            FILE *fL = fopen(filenameDensity, "rt");
                              for (i=0;i<par->numNodes; i++){
                                ierr = fscanf(fL, "%le %le",&par->nodos[i].nodalDensityInitial, &par->nodos[i].nodalDensity);
                              }
                            fclose(fL);
                            }

                        char filenamePressure[256];
                        sprintf(filenamePressure,"RestartPressure%d.dat",StepRestart);
                              if (  (fopen(filenamePressure, "rt")) == NULL)  {
                              PetscPrintf(PETSC_COMM_WORLD,"Pressure file NOT found!\n");
                              } else {
                                FILE *fL = fopen(filenamePressure, "rt");
                                for (i=0;i<par->numNodes; i++){
                                      ierr = fscanf(fL, "%le ",&par->nodos[i].nodalPressure);
                                }
                                fclose(fL);
                                }

  	  // Geometry

  		char filenameGeo[256];
  		sprintf(filenameGeo,"RestartGeo%d.dat",StepRestart);
  	  	    if (  (fopen(filenameGeo, "rt")) == NULL)  {
  	  	    PetscPrintf(PETSC_COMM_WORLD,"RestartGeo file NOT found!\n");
  	  	    } else {

  	  	    FILE *fL = fopen(filenameGeo, "rt");
            user->totalCurrentExplosiveVolume = 0;
            user->totalInitialExplosiveVolume = 0;
  	          for (i=0;i<par->numNodes; i++){

  	        	  ierr = fscanf(fL, "%le %le %le %le",&par->puntos[i].currentCoord[0], &par->puntos[i].currentCoord[1],&par->nodos[i].nodalVolumeInitial, &par->nodos[i].nodalVolume);
                 user->totalCurrentExplosiveVolume += par->nodos[i].nodalVolume;
                 user->totalInitialExplosiveVolume += par->nodos[i].nodalVolumeInitial;

  	          }
              PetscPrintf(PETSC_COMM_WORLD,"%e \n", user->totalCurrentExplosiveVolume);
              PetscPrintf(PETSC_COMM_WORLD,"%e \n", user->totalInitialExplosiveVolume);
  	  	    fclose(fL);
  	  	    }
  	  	    // GeometryOld
  	  	  PetscInt temp1 = StepRestart-1;
  	  		char filenameGeoOld[256];
  	  		sprintf(filenameGeoOld,"RestartGeo%d.dat",temp1);

  	  	  	    if (  (fopen(filenameGeoOld, "rt")) == NULL)  {
  	  	  	    PetscPrintf(PETSC_COMM_WORLD,"RestartGeoOld file NOT found!\n");
  	  	  	    } else {

  	  	  	    FILE *fL = fopen(filenameGeoOld, "rt");

  	  	          for (i=0;i<par->numNodes; i++){
  	  	        	  ierr = fscanf(fL, "%le %le",&par->puntos[i].totalPhysicalDisplacement[0], &par->puntos[i].totalPhysicalDisplacement[1]);
  	  	          }

  	  	  	    fclose(fL);
  	  	  	    }
  	  	  	for (i=0;i<par->numNodes; i++){
  	  	  	     par->puntos[i].totalPhysicalDisplacement[0] = par->puntos[i].currentCoord[0] - par->puntos[i].totalPhysicalDisplacement[0];
  	  	  	     par->puntos[i].totalPhysicalDisplacement[1] = par->puntos[i].currentCoord[1] - par->puntos[i].totalPhysicalDisplacement[1];
  	  	  	}
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "OutputOldGeometry"
PetscErrorCode OutputOldGeometry(PARAMETERS *par)
{

  PetscFunctionBegin;
  PetscErrorCode ierr;
  PetscInt i,j;
  PetscInt count=par->numNodes;
	POINTS *point = par->puntos;
  FILE *fileResult;
  char filenameResults[50];

     // ##################################################
     //                 Old Geometry File
     // ##################################################
     sprintf(filenameResults,"./RestartGeo%d.dat",par->stepNumber);
      fileResult=fopen(filenameResults,"wt");
      if (fileResult == NULL)
          {
      	printf("Error opening file!\n");
      	exit(1);
          }
        for(i=0; i<count;i++)
        {
         for(j=0; j<2;j++){
          PetscFPrintf(PETSC_COMM_WORLD,fileResult,"%.8e ", point[i].currentCoord[j]);
          if ((j+1) % 2 == 0) PetscFPrintf(PETSC_COMM_WORLD,fileResult,"\n  ");
        }
        }
      fclose(fileResult);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "OutputRestarts"
PetscErrorCode OutputRestarts(PARAMETERS *par,Vec U,Vec V)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  PetscInt i, j;

  PetscInt count=par->numNodes;
    MPI_Comm comm;
    PetscViewer viewer;
     ierr = PetscObjectGetComm((PetscObject)U,&comm);CHKERRQ(ierr);
     char filenameRestart[256];
     sprintf(filenameRestart,"RestartU%d.dat",par->stepNumber);
     ierr = PetscViewerBinaryOpen(comm,filenameRestart,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
     ierr = VecView(U,viewer);CHKERRQ(ierr);
     ierr = PetscObjectGetComm((PetscObject)V,&comm);CHKERRQ(ierr);
     sprintf(filenameRestart,"RestartV%d.dat",par->stepNumber);
     ierr = PetscViewerBinaryOpen(comm,filenameRestart,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
     ierr = VecView(V,viewer);CHKERRQ(ierr);
     ierr = PetscViewerDestroy(&viewer);
     FILE *fileResult;
     char filenameResults[50];

     // ##################################################
     //             Write Restart Time File
     // ##################################################
     sprintf(filenameResults,"./RestartTime%d.dat",par->stepNumber);
     fileResult=fopen(filenameResults,"wt");
     if (fileResult == NULL)
         {
     	printf("Error opening file!\n");
     	exit(1);
         }
         PetscFPrintf(PETSC_COMM_WORLD,fileResult,"%.8e  ", par->currentTime);

     fclose(fileResult);

  // ##################################################
  //                  Velocity File
  // ##################################################
     sprintf(filenameResults,"./RestartVelocity%d.dat",par->stepNumber);
     fileResult=fopen(filenameResults,"wt");
     if (fileResult == NULL)
         {
     	printf("Error opening file!\n");
     	exit(1);
         }
       for(i=0; i<count;i++)
       {
        for(j=0; j<2;j++){
         PetscFPrintf(PETSC_COMM_WORLD,fileResult,"%.8e  ", par->puntos[i].totalPhysicalVelocity[j]);
         if ((j+1) % 2 == 0) PetscFPrintf(PETSC_COMM_WORLD,fileResult,"\n  ");
       }
       }
     fclose(fileResult);

  // ##################################################
  //                  Acceleration File
  // ##################################################
      sprintf(filenameResults,"./RestartAcceleration%d.dat",par->stepNumber);
        fileResult=fopen(filenameResults,"wt");
        if (fileResult == NULL)
            {
        	printf("Error opening file!\n");
        	exit(1);
            }
          for(i=0; i<count;i++)
          {
           for(j=0; j<2;j++){
            PetscFPrintf(PETSC_COMM_WORLD,fileResult,"%.8e  ", par->puntos[i].totalPhysicalAcceleration[j]);
            if ((j+1) % 2== 0) PetscFPrintf(PETSC_COMM_WORLD,fileResult,"\n  ");
          }
          }
        fclose(fileResult);

     // ##################################################
     //                  Density File
     // ##################################################
    sprintf(filenameResults,"./RestartDensity%d.dat",par->stepNumber);
     fileResult=fopen(filenameResults,"wt");
     if (fileResult == NULL)
         {
     	printf("Error opening file!\n");
     	exit(1);
    }
       for(i=0; i<count;i++)
       {
         PetscFPrintf(PETSC_COMM_WORLD,fileResult,"%.8e  ", par->nodos[i].nodalDensityInitial);
         PetscFPrintf(PETSC_COMM_WORLD,fileResult,"%.8e  ", par->nodos[i].nodalDensity);
         PetscFPrintf(PETSC_COMM_WORLD,fileResult,"\n  ");
       }
     fclose(fileResult);

     // ##################################################
     //                  Pressure File
     // ##################################################
      sprintf(filenameResults,"./RestartPressure%d.dat",par->stepNumber);
       fileResult=fopen(filenameResults,"wt");
       if (fileResult == NULL)
           {
       	printf("Error opening file!\n");
       	exit(1);
           }
         for(i=0; i<count;i++)
         {
           PetscFPrintf(PETSC_COMM_WORLD,fileResult,"%.8e  \n", par->nodos[i].nodalPressure);
         }
       fclose(fileResult);

     // ##################################################
     //        Geometry & Associated Volume File
     // ##################################################
     sprintf(filenameResults,"./RestartGeo%d.dat",par->stepNumber);
      fileResult=fopen(filenameResults,"wt");
      if (fileResult == NULL)
          {
      	printf("Error opening file!\n");
      	exit(1);
          }
        for(i=0; i<count;i++)
        {
         for(j=0; j<2;j++){
          PetscFPrintf(PETSC_COMM_WORLD,fileResult,"%.8e  ", par->puntos[i].currentCoord[j]);
        }
          PetscFPrintf(PETSC_COMM_WORLD,fileResult,"%.8e  ", par->nodos[i].nodalVolumeInitial);
          PetscFPrintf(PETSC_COMM_WORLD,fileResult,"%.8e  \n", par->nodos[i].nodalVolume);
        }
      fclose(fileResult);
  PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "ComputeCurrentExplosiveVolume_dens"
PetscErrorCode ComputeCurrentExplosiveVolume_dens(AppCtx *user, PARAMETERS *par, const PetscScalar *arrayV,const PetscScalar *arrayU, PetscInt it)
{
// Update explosive volume and explosive nodal volume based on JWL EOS
//PetscPrintf(PETSC_COMM_WORLD,"  %e \n", u[0]);//Debug


  PetscReal density_initial = user->rhoRDX;
  PetscReal v0, v1;
	PetscErrorCode ierr;

  IGAElement         ele1;
  ierr = IGAElementCreate(&ele1);CHKERRQ(ierr);
  ierr = IGAElementInit(ele1,user->iga);CHKERRQ(ierr);

  PetscScalar *V;
  PetscScalar *U;
  // Get the pressures to output at the same time:
  PetscReal dens0  = 1659.0;
  PetscReal P0     = 100000.0;
  PetscReal A      = 495.1e9;
  PetscReal B      = 7.21e9;
  PetscReal C      = 1.62e9;
  PetscReal R1     = 4.387;
  PetscReal R2     = 0.9954;
  PetscReal omega  = 0.3469;
  PetscReal E0     = 5877.9e3;
  PetscReal Pcr    = 2.0e11;
  PetscReal nu     = 0;
  PetscReal Ptest  = A*exp(-R1*nu) + B*exp(-R2*nu) + C/(pow(nu, 1+omega));

  PetscReal P        = Pcr;
  PetscReal fprime   = 0.0;
  PetscReal cs       = 0.0;
  user->totalCurrentExplosiveVolume = 0.0;

  PetscReal q[2];
	   for (int i=0; i<(par->numNodes); i++){
          ele1->nval = 0;
    			q[0] = par->puntos[i].currentCoord[0]/user->Lx;
    			q[1] = par->puntos[i].currentCoord[1]/user->Ly;

          IGAPoint pnt1;
          ierr = IGAPointCreate(&pnt1);CHKERRQ(ierr);
            pnt1->dof   = user->iga->dof;
    				pnt1->dim   = user->iga->dim;
    				pnt1->point = q;

	  if (par->puntos[i].material == 1){
    if (IGALocateElement(user->iga,pnt1->point,ele1)) {
      ierr = IGAElementGetWorkVal(ele1,&U);CHKERRQ(ierr);
      ierr = IGAElementGetValues(ele1,arrayU,U);CHKERRQ(ierr);
      ierr = IGAPointInit(pnt1,ele1);CHKERRQ(ierr);
      ierr = IGAPointEval(user->iga,pnt1);CHKERRQ(ierr);

    	    	  PetscScalar u[pnt1->dof];
    	    	  IGAPointFormValue(pnt1,U,&u[0]);

              PetscScalar dens = u[0];
                if(par->stepNumber == 0){
                  par->nodos[i].nodalDensityInitial = dens;
                }
                par->nodos[i].nodalDensity = dens;
                PetscScalar ux  = u[1];
                PetscScalar uy  = u[2];
                PetscScalar temp= u[3];
              if(dens<=0.0){
                PetscPrintf(PETSC_COMM_WORLD,"Density < 0, Current Explosive Volume Update Step error \n");
                exit(0);}

              if(par->stepNumber>0){
                par->nodos[i].nodalVolume = par->nodos[i].nodalVolumeInitial*par->nodos[i].nodalDensityInitial/dens;
                user->totalCurrentExplosiveVolume += par->nodos[i].nodalVolume;
              }


              if(dens<0){exit(0);}
              nu     =  par->nodos[i].nodalDensity/dens0;
              Ptest  = A*exp(-R1*nu) + B*exp(-R2*nu) + C/(pow(nu, 1+omega));
              P      = Pcr;

              if(Ptest > Pcr){
                P        = Ptest;
              }else{
                P        = A*(1-omega/(R1*nu))*exp(-R1*nu) + B*(1-omega/(R2*nu))*exp(-R2*nu) + omega*dens0*E0/nu;
              }
              par->nodos[i].nodalPressure = P;
              }
        ierr = IGAPointDestroy(&pnt1);CHKERRQ(ierr);
				  }
	  	  }
    ierr = IGAElementDestroy(&ele1);CHKERRQ(ierr);

	  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeDeformationGradient"
PetscErrorCode ComputeDeformationGradient(PARAMETERS *par,PetscReal af)
{
	  PetscErrorCode ierr;

	  for(int i = 0; i<par->numNodes; i++){

			  PetscReal temp1[4] = {0.0};
			  PetscReal temp2[4] = {0.0};
			  PetscReal temp3;
			  PetscReal inv[4]	 = {0.0};
			  PetscReal inter[4] = {0.0};
			  PetscReal	temp = 0.0;

			  temp1[0] = 1.0 - par->timeStep*par->nodos[i].velocityGradient[0]*af;
			  temp1[1] = -par->timeStep*par->nodos[i].velocityGradient[1]*af;
			  temp1[2] = -par->timeStep*par->nodos[i].velocityGradient[2]*af;
			  temp1[3] = 1.0 -par->timeStep*par->nodos[i].velocityGradient[3]*af;

			  temp2[0] = 1.0 + par->timeStep*par->nodos[i].velocityGradient[0]*(1-af);
			  temp2[1] = par->timeStep*par->nodos[i].velocityGradient[1]*(1-af);
			  temp2[2] = par->timeStep*par->nodos[i].velocityGradient[2]*(1-af);
			  temp2[3] = 1.0 + par->timeStep*par->nodos[i].velocityGradient[3]*(1-af);

        temp = temp1[0]*temp1[3] - temp1[1]*temp1[2];

			  inv[0] = temp1[3]/temp;
			  inv[1] = -temp1[1]/temp;
			  inv[2] = -temp1[2]/temp;
			  inv[3] = temp1[0]/temp;

			  inter[0] = inv[0]*temp2[0] + inv[1]*temp2[2];
			  inter[1] = inv[0]*temp2[1] + inv[1]*temp2[3];
			  inter[2] = inv[2]*temp2[0] + inv[3]*temp2[2];
			  inter[3] = inv[2]*temp2[1] + inv[3]*temp2[3];
        //PetscPrintf(PETSC_COMM_WORLD,"%e %e \n", inter[0], inter[1]);
			  temp1[0] = inter[0]*par->nodos[i].DeformationGradientOld[0] + inter[1]*par->nodos[i].DeformationGradientOld[2];
			  temp1[1] = inter[0]*par->nodos[i].DeformationGradientOld[1] + inter[1]*par->nodos[i].DeformationGradientOld[3];
			  temp1[2] = inter[2]*par->nodos[i].DeformationGradientOld[0] + inter[3]*par->nodos[i].DeformationGradientOld[2];
			  temp1[3] = inter[2]*par->nodos[i].DeformationGradientOld[1] + inter[3]*par->nodos[i].DeformationGradientOld[3];

			  par->nodos[i].alphaDeformationGradient[0] = par->nodos[i].DeformationGradientOld[0] + af*(temp1[0] - par->nodos[i].DeformationGradientOld[0]);
			  par->nodos[i].alphaDeformationGradient[1] = par->nodos[i].DeformationGradientOld[1] + af*(temp1[1] - par->nodos[i].DeformationGradientOld[1]);
			  par->nodos[i].alphaDeformationGradient[2] = par->nodos[i].DeformationGradientOld[2] + af*(temp1[2] - par->nodos[i].DeformationGradientOld[2]);
			  par->nodos[i].alphaDeformationGradient[3] = par->nodos[i].DeformationGradientOld[3] + af*(temp1[3] - par->nodos[i].DeformationGradientOld[3]);

			  par->nodos[i].currentDeformationGradient[0] = temp1[0];
			  par->nodos[i].currentDeformationGradient[1] = temp1[1];
			  par->nodos[i].currentDeformationGradient[2] = temp1[2];
			  par->nodos[i].currentDeformationGradient[3] = temp1[3];

			  par->nodos[i].determinantAlphaDeformationGradient = par->nodos[i].alphaDeformationGradient[0]*par->nodos[i].alphaDeformationGradient[3]-par->nodos[i].alphaDeformationGradient[1]*par->nodos[i].alphaDeformationGradient[2];
	 	}

	  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ComputeVelocityGradient"
PetscErrorCode ComputeVelocityGradient(PARAMETERS *par, AppCtx *user,Vec vecU,const PetscScalar *arrayU)

{

  PetscErrorCode ierr;
  PetscScalar       *U;
  IGAElement ele1;
  IGAElementCreate(&ele1);
  IGAElementInit(ele1,user->iga);

   for(int i = 0; i<par->numNodes; i++){
				  IGAPoint           pnt1;
				  ierr = IGAPointCreate(&pnt1);CHKERRQ(ierr);
						PetscReal q1[2];
						q1[0] = par->puntos[i].currentCoord[0]/user->Lx;
						q1[1] = par->puntos[i].currentCoord[1]/user->Ly;
						ierr = IGAPointReset(pnt1);CHKERRQ(ierr);
						pnt1->dof   = user->iga->dof;
            pnt1->dim   = user->iga->dim;
						pnt1->point = q1;
	          ele1->nval = 0;
				if (IGALocateElement(user->iga,pnt1->point,ele1)) {
								ierr = IGAPointInit(pnt1,ele1);CHKERRQ(ierr);
								ierr = IGAPointEval(user->iga,pnt1);CHKERRQ(ierr);
								ierr = IGAElementGetWorkVal(ele1,&U);CHKERRQ(ierr);
								ierr = IGAElementGetValues(ele1,arrayU,U);CHKERRQ(ierr);
								PetscScalar grad_u[4][2];
								IGAPointFormGrad(pnt1,U,&grad_u[0][0]);
								par->nodos[i].velocityGradient[0] = grad_u[1][0];
								par->nodos[i].velocityGradient[1] = grad_u[1][1];
								par->nodos[i].velocityGradient[2] = grad_u[2][0];
								par->nodos[i].velocityGradient[3] = grad_u[2][1];
				}
        ierr = IGAPointDestroy(&pnt1);CHKERRQ(ierr);
		   }
  ierr = IGAElementDestroy(&ele1);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef  __FUNCT__
#define __FUNCT__ "IGAComputeResidualFS"
PetscErrorCode IGAComputeResidualFS(PARAMETERS *par,IGA iga,PetscReal dt,
        PetscReal a,const PetscScalar *arrayV,
        PetscReal t,const PetscScalar *arrayU,
        Vec vecF,AppCtx *user)
{

  PetscScalar       *V,*U,*Kvec,*Kvec1;
  PetscErrorCode    ierr;
  IGAElement        ele;

  /* Clear global matrix J*/
   ierr = VecZeroEntries(vecF);CHKERRQ(ierr);
   ierr = IGAElementCreate(&ele);CHKERRQ(ierr);
   ierr = IGAElementInit(ele,user->iga);CHKERRQ(ierr);

   PetscInt i;
   PetscInt j;
   PetscReal q[2];
for(j = 0; j<par->numNodes; j++){
    ele->nval = 0;
 		q[0] = par->puntos[j].currentCoord[0]/user->Lx;
 		q[1] = par->puntos[j].currentCoord[1]/user->Ly;

    IGAPoint pnt;
    ierr = IGAPointCreate(&pnt);CHKERRQ(ierr);
    pnt->dof   = user->iga->dof;
    pnt->dim   = user->iga->dim;
    pnt->point = q;

 			if (IGALocateElement(user->iga,pnt->point,ele)) {

        ierr = IGAPointInit(pnt,ele);CHKERRQ(ierr);
        ierr = IGAPointEval(user->iga,pnt);CHKERRQ(ierr);

 			 	ierr = IGAElementGetWorkVal(ele,&V);CHKERRQ(ierr);
 			 	ierr = IGAElementGetWorkVal(ele,&U);CHKERRQ(ierr);
 			 	ierr = IGAElementGetValues(ele,arrayV,V);CHKERRQ(ierr);
 			 	ierr = IGAElementGetValues(ele,arrayU,U);CHKERRQ(ierr);

 	    	ierr = IGAPointGetWorkVec(pnt,&Kvec);CHKERRQ(ierr);
 	    	ierr = IGAPointGetWorkVec(pnt,&Kvec1);CHKERRQ(ierr);

 	    	ierr = ResidualFS(pnt,dt,a,V,t,U,Kvec,user);CHKERRQ(ierr);
 	    	ierr = ResidualRDX(pnt,dt,a,V,t,U,Kvec1,user);CHKERRQ(ierr);

 	    	for(i=0; i< (ele->nen)*(ele->dof); i++){
 	    	        	 Kvec[i] = (Kvec1[i] - Kvec[i])*par->nodos[j].nodalVolume;
 	    	}
 	    	ierr = IGAElementAssembleVec(ele,Kvec,vecF);CHKERRQ(ierr);
        ierr = IGAPointDestroy(&pnt);CHKERRQ(ierr);
 			}
     }
  ierr = IGAElementDestroy(&ele);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "interpolateVelocityOnSolidNodes"
PetscErrorCode interpolateVelocityOnSolidNodes(PARAMETERS *par, AppCtx *user,PetscInt in, Vec vecU,const PetscScalar *arrayU)
{

  PetscErrorCode     ierr;
  PetscScalar        *U;
  PetscScalar        uf[4];
  PetscReal          pt[2];
  IGAElement         eleF;
  IGAPoint           pntF;

  ierr = IGAPointCreate(&pntF);CHKERRQ(ierr);
  ierr = IGAElementCreate(&eleF);CHKERRQ(ierr);
  ierr = IGAElementInit(eleF,user->iga);CHKERRQ(ierr);

  pt[0] =  par->puntos[in].currentCoord[0]/user->Lx;
  pt[1] =  par->puntos[in].currentCoord[1]/user->Ly;
	pntF->dof   = user->iga->dof;
	pntF->dim   = user->iga->dim;
	pntF->point = pt;

		if (IGALocateElement(user->iga,pntF->point,eleF)) {

		  ierr = IGAElementGetWorkVal(eleF,&U);CHKERRQ(ierr);
	  	ierr = IGAElementGetValues(eleF,arrayU,U);CHKERRQ(ierr);
		  ierr = IGAPointInit(pntF,eleF);CHKERRQ(ierr);
  		ierr = IGAPointEval(user->iga,pntF);CHKERRQ(ierr);
		  IGAPointFormValue(pntF,U,&uf[0]);

		    par->puntos[in].totalPhysicalVelocity[0] = uf[1];
		    par->puntos[in].totalPhysicalVelocity[1] = uf[2];
		}

	    ierr = IGAElementDestroy(&eleF);CHKERRQ(ierr);
	    ierr = IGAPointDestroy(&pntF);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


#undef  __FUNCT__
#define __FUNCT__ "IGAComputeIJacobianFS"
PetscErrorCode IGAComputeIJacobianFS(PetscReal shift,Mat MassFS,AppCtx *user,const PetscScalar *arrayU, PARAMETERS *par)
{


  PetscScalar       *U;
  PetscErrorCode    ierr;
  IGAElement         ele;

  PetscReal        A0[4][4]       = {{0.0}};
  PetscReal        A01[4][4]      = {{0.0}};
  PetscScalar      KK[16];
  PetscScalar      KK1[16];
  ierr = IGAElementCreate(&ele);CHKERRQ(ierr);
  ierr = IGAElementInit(ele,user->iga);CHKERRQ(ierr);

  PetscInt i;
  PetscReal q[2];
    for(i = 0; i<par->numNodes; i++){
      ele->nval = 0;
			q[0] = par->puntos[i].currentCoord[0]/user->Lx;
			q[1] = par->puntos[i].currentCoord[1]/user->Ly;

      IGAPoint pnt;
      ierr = IGAPointCreate(&pnt);CHKERRQ(ierr);
        pnt->dof   = user->iga->dof;
				pnt->dim   = user->iga->dim;
				pnt->point = q;

		   if (IGALocateElement(user->iga,pnt->point,ele)) {
         ierr = IGAElementGetWorkVal(ele,&U);CHKERRQ(ierr);
         ierr = IGAElementGetValues(ele,arrayU,U);CHKERRQ(ierr);
         ierr = IGAPointInit(pnt,ele);CHKERRQ(ierr);
         ierr = IGAPointEval(user->iga,pnt);CHKERRQ(ierr);


				PetscReal *N0 = pnt->shape[0];
				PetscInt  nen = pnt->nen;
				PetscInt  dof = pnt->dof;
				PetscInt  a;
        PetscReal Cv     = user->Cv;

				PetscScalar u[dof];
				IGAPointFormValue(pnt,U,&u[0]);

        PetscScalar dens= u[0];
        PetscScalar ux  = u[1];
        PetscScalar uy  = u[2];
        PetscScalar temp= u[3];

        A0[0][0] = 1.0;
        A0[1][0] = ux;
        A0[1][1] = dens;
        A0[2][0] = uy;
        A0[2][2] = dens;
        A0[3][0] = Cv*temp;
        A0[3][3] = dens*Cv;

        Cv = 143.3; //this is not the right Cv for RDX in this situation
        //Just use this for now...

        A01[0][0] = 1.0;
        A01[1][0] = ux;
        A01[1][1] = dens;
        A01[2][0] = uy;
        A01[2][2] = dens;
        A01[3][0] = Cv*temp;
        A01[3][3] = dens*Cv;

				for (a=0; a<nen; a++) {

					KK[0] = A0[0][0]*shift*N0[a]*par->nodos[i].nodalVolume;
					KK[1] = A0[0][1]*shift*N0[a]*par->nodos[i].nodalVolume;
					KK[2] = A0[0][2]*shift*N0[a]*par->nodos[i].nodalVolume;
					KK[3] = A0[0][3]*shift*N0[a]*par->nodos[i].nodalVolume;

					KK[4] = A0[1][0]*shift*N0[a]*par->nodos[i].nodalVolume;
					KK[5] = A0[1][1]*shift*N0[a]*par->nodos[i].nodalVolume;
					KK[6] = A0[1][2]*shift*N0[a]*par->nodos[i].nodalVolume;
					KK[7] = A0[1][3]*shift*N0[a]*par->nodos[i].nodalVolume;

					KK[8] = A0[2][0]*shift*N0[a]*par->nodos[i].nodalVolume;
					KK[9] = A0[2][1]*shift*N0[a]*par->nodos[i].nodalVolume;
					KK[10] = A0[2][2]*shift*N0[a]*par->nodos[i].nodalVolume;
					KK[11] = A0[2][3]*shift*N0[a]*par->nodos[i].nodalVolume;

					KK[12] = A0[3][0]*shift*N0[a]*par->nodos[i].nodalVolume;
					KK[13] = A0[3][1]*shift*N0[a]*par->nodos[i].nodalVolume;
					KK[14] = A0[3][2]*shift*N0[a]*par->nodos[i].nodalVolume;
					KK[15] = A0[3][3]*shift*N0[a]*par->nodos[i].nodalVolume;

					KK1[0] = A01[0][0]*shift*N0[a]*par->nodos[i].nodalVolume;
					KK1[1] = A01[0][1]*shift*N0[a]*par->nodos[i].nodalVolume;
					KK1[2] = A01[0][2]*shift*N0[a]*par->nodos[i].nodalVolume;
					KK1[3] = A01[0][3]*shift*N0[a]*par->nodos[i].nodalVolume;

					KK1[4] = A01[1][0]*shift*N0[a]*par->nodos[i].nodalVolume;
					KK1[5] = A01[1][1]*shift*N0[a]*par->nodos[i].nodalVolume;
					KK1[6] = A01[1][2]*shift*N0[a]*par->nodos[i].nodalVolume;
					KK1[7] = A01[1][3]*shift*N0[a]*par->nodos[i].nodalVolume;

					KK1[8] = A01[2][0]*shift*N0[a]*par->nodos[i].nodalVolume;
					KK1[9] = A01[2][1]*shift*N0[a]*par->nodos[i].nodalVolume;
					KK1[10] = A01[2][2]*shift*N0[a]*par->nodos[i].nodalVolume;
					KK1[11] = A01[2][3]*shift*N0[a]*par->nodos[i].nodalVolume;

					KK1[12] = A01[3][0]*shift*N0[a]*par->nodos[i].nodalVolume;
					KK1[13] = A01[3][1]*shift*N0[a]*par->nodos[i].nodalVolume;
					KK1[14] = A01[3][2]*shift*N0[a]*par->nodos[i].nodalVolume;
					KK1[15] = A01[3][3]*shift*N0[a]*par->nodos[i].nodalVolume;


					KK[0] -= KK1[0];
					KK[1] -= KK1[1];
					KK[2] -= KK1[2];
					KK[3] -= KK1[3];
					KK[4] -= KK1[4];
					KK[5] -= KK1[5];
					KK[6] -= KK1[6];
					KK[7] -= KK1[7];
					KK[8] -= KK1[8];
					KK[9] -= KK1[9];
					KK[10] -= KK1[10];
					KK[11] -= KK1[11];
					KK[12] -= KK1[12];
					KK[13] -= KK1[13];
					KK[14] -= KK1[14];
					KK[15] -= KK1[15];

					PetscInt ID = ele->mapping[a];
					PetscInt index_array[4] = {0};
					index_array[0] = ID*dof;
					index_array[1] = ID*dof + 1;
					index_array[2] = ID*dof + 2;
					index_array[3] = ID*dof + 3;
					MatSetValuesLocal(MassFS,4,index_array,4,index_array,KK,ADD_VALUES);CHKERRQ(ierr);
					}
          ierr = IGAPointDestroy(&pnt);CHKERRQ(ierr);
		   }
    }
  ierr = MatAssemblyBegin(MassFS,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd  (MassFS,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = IGAElementDestroy(&ele);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ComputeCurrentExplosiveVolume"
PetscErrorCode ComputeCurrentExplosiveVolume(AppCtx *user, PARAMETERS *par, const PetscScalar *arrayV,const PetscScalar *arrayU, PetscInt it)
{
PetscErrorCode ierr;
  user->totalCurrentExplosiveVolume = 0.0;
	  	  for (int i=0; i<(par->numNodes); i++){
	  		  if (par->puntos[i].material == 1){
            if(par->stepNumber==0){
              par->nodos[i].nodalVolumeInitial = par->nodos[i].nodalVolume;
            }
            //PetscPrintf(PETSC_COMM_WORLD,"%e \n", par->nodos[i].determinantAlphaDeformationGradient);
            par->nodos[i].nodalVolume = par->nodos[i].nodalVolumeInitial*par->nodos[i].determinantAlphaDeformationGradient;
            user->totalCurrentExplosiveVolume += par->nodos[i].nodalVolume;
          }
       }
       PetscReal density_initial = user->rhoRDX;
       PetscReal v0;


       //user->totalCurrentExplosiveVolume = 0.0;
       IGAElement         ele1;
       ierr = IGAElementCreate(&ele1);CHKERRQ(ierr);
       IGAPoint           pnt1;
       ierr = IGAPointCreate(&pnt1);CHKERRQ(ierr);
       PetscScalar *V;
       PetscScalar *U;
       // Get the pressures to output at the same time:
       PetscReal dens0  = 1659.0;
       PetscReal P0     = 100000.0;
       PetscReal A      = 495.1e9;
       PetscReal B      = 7.21e9;
       PetscReal C      = 1.62e9;
       PetscReal R1     = 4.387;
       PetscReal R2     = 0.9954;
       PetscReal omega  = 0.3469;
       PetscReal E0     = 5877.9e3;
       PetscReal Pcr    = 2.0e11;

       PetscReal nu     = 0;
       PetscReal Ptest  = A*exp(-R1*nu) + B*exp(-R2*nu) + C/(pow(nu, 1+omega));
       PetscReal P        = Pcr;
       PetscReal fprime   = 0.0;
       PetscReal cs       = 0.0;

             for (int i=0; i<(par->numNodes); i++){

               POINTS *pointSolid = &par->puntos[i];
               PetscReal q1[2];
               q1[0] = par->puntos[i].currentCoord[0]/user->Lx;
               q1[1] = par->puntos[i].currentCoord[1]/user->Ly;
               ierr = IGAPointReset(pnt1);CHKERRQ(ierr);

               if (par->puntos[i].material == 1){
                 ierr = IGAPointCreate(&pnt1);CHKERRQ(ierr);
                 ierr = IGAElementCreate(&ele1);CHKERRQ(ierr);

                 ierr = IGAPointReset(pnt1);CHKERRQ(ierr);
                 ierr = IGAElementInit(ele1,user->iga);CHKERRQ(ierr);
                 pnt1->dof   = user->iga->dof;
                 pnt1->dim   = user->iga->dim;
                 pnt1->point = q1;

                 if (IGALocateElement(user->iga,pnt1->point,ele1)) {

                   ierr = IGAPointInit(pnt1,ele1);CHKERRQ(ierr);
                   ierr = IGAPointEval(user->iga,pnt1);CHKERRQ(ierr);
                   PetscInt nen = ele1->nen;
                   PetscInt dim = ele1->dim;

                   ierr = IGAElementGetWorkVal(ele1,&V);CHKERRQ(ierr);
                   ierr = IGAElementGetValues(ele1,arrayV,V);CHKERRQ(ierr);
                   PetscScalar v[pnt1->dof];
                   IGAPointFormValue(pnt1,V,&v[0]);

                   ierr = IGAElementGetWorkVal(ele1,&U);CHKERRQ(ierr);
                   ierr = IGAElementGetValues(ele1,arrayU,U);CHKERRQ(ierr);
                     PetscScalar grad_u[pnt1->dof][pnt1->dim];
                     IGAPointFormGrad (pnt1,U,&grad_u[0][0]);
                     PetscScalar u[pnt1->dof];
                     IGAPointFormValue(pnt1,U,&u[0]);

                   //par->nodos[i].nodalVolume = u[0]*par->nodos[i].nodalVolumeInitial/density_initial;
                   //user->totalCurrentExplosiveVolume += par->nodos[i].nodalVolume;
                   par->nodos[i].nodalDensity = u[0];
                   nu     =  par->nodos[i].nodalDensity/dens0;
                   Ptest  = A*exp(-R1*nu) + B*exp(-R2*nu) + C/(pow(nu, 1+omega));
                   P      = Pcr;
                   if(Ptest > Pcr){
                     P        = Ptest;
                   }else{
                     P        = A*(1-omega/(R1*nu))*exp(-R1*nu) + B*(1-omega/(R2*nu))*exp(-R2*nu) + omega*dens0*E0/nu;
                   }
                   par->nodos[i].nodalPressure = P;
                   ierr = IGAElementDestroy(&ele1);CHKERRQ(ierr);
                   ierr = IGAPointDestroy(&pnt1);CHKERRQ(ierr);
                       }
               }
             }

             ierr = IGAElementDestroy(&ele1);CHKERRQ(ierr);
             ierr = IGAPointDestroy(&pnt1);CHKERRQ(ierr);

	  	PetscFunctionReturn(0);
}


#undef  __FUNCT__
#define __FUNCT__ "input"
//This function reads the input files for the particle data, specifically foreground.dat **2d specific**
// Read file structured as:
// #nodes 0 0 0 0 0 0 0
// ID x y z dx dy dz vol
PetscErrorCode input (PARAMETERS *par, AppCtx *user)
{
	PetscErrorCode  ierr;
	PetscInt i,j;
  PetscReal temp = 0;
  PetscPrintf(PETSC_COMM_WORLD,"Reading Input\n");
	//--------------------
    //  RDX.dat
	//--------------------
	//Here we read the foreground.dat file and we assign the data of the file to the initialCoord of the particles.
	//Obviously if the computation starts at time zero these coordinates are also the currentCoord

	char filenameCoord[256];
	sprintf(filenameCoord,"./RDX.dat");
    	    if ((fopen(filenameCoord, "rt")) == NULL){
    	    PetscPrintf(PETSC_COMM_WORLD,"RDX.dat file NOT found!\n");
    	    } else {
    	    FILE *fL = fopen(filenameCoord, "rt");
    	    ierr = fscanf(fL, "%d %le %le %le %le %le %le %le", &par->numNodes, &temp, &temp, &temp, &temp, &temp, &temp, &temp);

    	    par->numPoints = par->numNodes;

    	    PetscPrintf(PETSC_COMM_WORLD," numNodes %d \n",par->numNodes);
    	    ierr = PetscMalloc1(sizeof(*par->puntos)*par->numPoints,POINTS,&par->puntos);CHKERRQ(ierr);
    	    ierr = PetscMalloc1(sizeof(*par->nodos)*par->numNodes,NODE,&par->nodos);CHKERRQ(ierr);

    	    //Here is where we read from the file with the command fscanf
                user->totalInitialExplosiveVolume=0;
                for (i=0; i<(par->numNodes); i++){
                  ierr = fscanf(fL, "%d %le %le %le %le %le %le %le",
                  &par->puntos[i].ID, &par->puntos[i].initialCoord[0],
                  &par->puntos[i].initialCoord[1], &temp, &par->puntos[i].hvect[0], &par->puntos[i].hvect[1],
                  &temp, &par->nodos[i].nodalVolume);

                par->puntos[i].ID=par->puntos[i].ID-1;
                par->nodos[i].nodalVolumeInitial = par->nodos[i].nodalVolume;

                user->totalInitialExplosiveVolume+=par->nodos[i].nodalVolume;
                user->totalCurrentExplosiveVolume+=par->nodos[i].nodalVolume;
                par->puntos[i].material = 1; //material 1 corresponds to RDX

                par->puntos[i].initialCoord[0] = par->puntos[i].initialCoord[0]+0.25/2.0-0.005/2.0;
                }
                fclose(fL);
    	    }

    	    for (i = 0; i < (par->numNodes); i++){
    	      par->nodos[i].nodeID = par->puntos[i].ID;
    	    	//Here we set the current coordinates to be the same as the initial coordinates, which is true only at the beggining of the computation
    	    	par->puntos[i].currentCoord[0] = par->puntos[i].initialCoord[0];
    	    	par->puntos[i].currentCoord[1] = par->puntos[i].initialCoord[1];

		       //PetscPrintf(PETSC_COMM_WORLD,"%le \n", par->puntos[par->numNodes].initialCoord[1]);
    	    //Here we are initializing the displacement and velocity of the particles at the beggining to be zero
            for (j=0;j<2; j++){
            par->puntos[i].totalPhysicalDisplacement[j] = 0.0;
            par->puntos[i].totalPhysicalVelocity[j] = 0.0;
            }
            //Same for the stresses and strains
            for (j=0;j<3; j++){
            par->puntos[i].totalStrain[j] = 0.0;
            par->puntos[i].totalStress[j] = 0.0;
            par->puntos[i].totalStressderivativeX[j]=0.0;
            par->puntos[i].totalStressderivativeY[j]=0.0;
            }
  }
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "Compute_inverse"
PetscErrorCode Compute_inverse(PetscReal (*matrix)[4],PetscReal (*matrix_inv)[4])
{

	PetscReal m[16] = {0.0};
	PetscReal invOut[16] = {0.0};
	PetscInt i,j;

	for (i=0;i<4;i++)
		for (j=0;j<4;j++)
			m[4*i+j] = matrix[i][j];

	PetscReal inv[16], det;


    inv[0] = m[5]  * m[10] * m[15] -
             m[5]  * m[11] * m[14] -
             m[9]  * m[6]  * m[15] +
             m[9]  * m[7]  * m[14] +
             m[13] * m[6]  * m[11] -
             m[13] * m[7]  * m[10];

    inv[4] = -m[4]  * m[10] * m[15] +
              m[4]  * m[11] * m[14] +
              m[8]  * m[6]  * m[15] -
              m[8]  * m[7]  * m[14] -
              m[12] * m[6]  * m[11] +
              m[12] * m[7]  * m[10];

    inv[8] = m[4]  * m[9] * m[15] -
             m[4]  * m[11] * m[13] -
             m[8]  * m[5] * m[15] +
             m[8]  * m[7] * m[13] +
             m[12] * m[5] * m[11] -
             m[12] * m[7] * m[9];

    inv[12] = -m[4]  * m[9] * m[14] +
               m[4]  * m[10] * m[13] +
               m[8]  * m[5] * m[14] -
               m[8]  * m[6] * m[13] -
               m[12] * m[5] * m[10] +
               m[12] * m[6] * m[9];

    inv[1] = -m[1]  * m[10] * m[15] +
              m[1]  * m[11] * m[14] +
              m[9]  * m[2] * m[15] -
              m[9]  * m[3] * m[14] -
              m[13] * m[2] * m[11] +
              m[13] * m[3] * m[10];

    inv[5] = m[0]  * m[10] * m[15] -
             m[0]  * m[11] * m[14] -
             m[8]  * m[2] * m[15] +
             m[8]  * m[3] * m[14] +
             m[12] * m[2] * m[11] -
             m[12] * m[3] * m[10];

    inv[9] = -m[0]  * m[9] * m[15] +
              m[0]  * m[11] * m[13] +
              m[8]  * m[1] * m[15] -
              m[8]  * m[3] * m[13] -
              m[12] * m[1] * m[11] +
              m[12] * m[3] * m[9];

    inv[13] = m[0]  * m[9] * m[14] -
              m[0]  * m[10] * m[13] -
              m[8]  * m[1] * m[14] +
              m[8]  * m[2] * m[13] +
              m[12] * m[1] * m[10] -
              m[12] * m[2] * m[9];

    inv[2] = m[1]  * m[6] * m[15] -
             m[1]  * m[7] * m[14] -
             m[5]  * m[2] * m[15] +
             m[5]  * m[3] * m[14] +
             m[13] * m[2] * m[7] -
             m[13] * m[3] * m[6];

    inv[6] = -m[0]  * m[6] * m[15] +
              m[0]  * m[7] * m[14] +
              m[4]  * m[2] * m[15] -
              m[4]  * m[3] * m[14] -
              m[12] * m[2] * m[7] +
              m[12] * m[3] * m[6];

    inv[10] = m[0]  * m[5] * m[15] -
              m[0]  * m[7] * m[13] -
              m[4]  * m[1] * m[15] +
              m[4]  * m[3] * m[13] +
              m[12] * m[1] * m[7] -
              m[12] * m[3] * m[5];

    inv[14] = -m[0]  * m[5] * m[14] +
               m[0]  * m[6] * m[13] +
               m[4]  * m[1] * m[14] -
               m[4]  * m[2] * m[13] -
               m[12] * m[1] * m[6] +
               m[12] * m[2] * m[5];

    inv[3] = -m[1] * m[6] * m[11] +
              m[1] * m[7] * m[10] +
              m[5] * m[2] * m[11] -
              m[5] * m[3] * m[10] -
              m[9] * m[2] * m[7] +
              m[9] * m[3] * m[6];

    inv[7] = m[0] * m[6] * m[11] -
             m[0] * m[7] * m[10] -
             m[4] * m[2] * m[11] +
             m[4] * m[3] * m[10] +
             m[8] * m[2] * m[7] -
             m[8] * m[3] * m[6];

    inv[11] = -m[0] * m[5] * m[11] +
               m[0] * m[7] * m[9] +
               m[4] * m[1] * m[11] -
               m[4] * m[3] * m[9] -
               m[8] * m[1] * m[7] +
               m[8] * m[3] * m[5];

    inv[15] = m[0] * m[5] * m[10] -
              m[0] * m[6] * m[9] -
              m[4] * m[1] * m[10] +
              m[4] * m[2] * m[9] +
              m[8] * m[1] * m[6] -
              m[8] * m[2] * m[5];

    det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];


    det = 1.0 / det;

    for (i = 0; i < 16; i++)
        invOut[i] = inv[i] * det;


	for (i=0;i<4;i++)
		for (j=0;j<4;j++)
			matrix_inv[i][j] = invOut[4*i+j];

return 0;
}

#undef __FUNCT__
#define __FUNCT__ "Compute_square_root_inverse"
PetscErrorCode Compute_square_root_inverse(PetscReal (*TauTemp)[4],PetscReal (*TauTempInv)[4],PetscReal taut,void *ctx)
{

    AppCtx *user = (AppCtx *)ctx;

//     PetscInt  dim    = user->iga->dim;
     PetscInt  dof    = user->iga->dof;

     PetscInt i,j,l,aa;


     PetscReal Y[4][4]={{0.0}};
     PetscReal Z[4][4]={{0.0}};
     PetscReal Yinv[4][4]={{0.0}};
     PetscReal Zinv[4][4]={{0.0}};

     for (i=0;i<dof;i++){

         Z[i][i] = 1.0/sqrt(taut);

         for (j=0;j<dof;j++){
             Y[i][j]= TauTemp[i][j];
         }
     }




     for (aa=0;aa<10;aa++){



    	 Compute_inverse(Y,Yinv);
    	 Compute_inverse(Z,Zinv);





          for (i=0;i<dof;i++){
              for (j=0;j<dof;j++){
                  Y[i][j]= 0.5*(Y[i][j] + Zinv[i][j]);
              }
          }

          for (i=0;i<dof;i++){
              for (j=0;j<dof;j++){
                  Z[i][j]= 0.5*(Z[i][j] + Yinv[i][j]);
              }
          }

     }


     for (i=0;i<dof;i++){
         for (j=0;j<dof;j++){
             TauTempInv[i][j]= Z[i][j];
         }
     }


  return 0;
}




#undef __FUNCT__
#define __FUNCT__ "DirectTau"
PetscErrorCode DirectTau(PetscReal G[2][2],
         PetscReal dt,PetscScalar u[],
         PetscScalar (*tau)[4],void *ctx,PetscReal (*A0inv)[4],PetscReal (*A1_cons)[4],PetscReal (*A2_cons)[4],PetscReal (*K_cons)[2][4][4],PetscReal *umi)
{

    AppCtx *user = (AppCtx *)ctx;

     PetscInt  dim    = user->iga->dim;

     PetscInt  dof    = user->iga->dof;

     PetscInt i,j;



  //Tau for Conservation variables


   PetscReal TauTemp[4][4]={{0.0}};
   PetscReal TauTempInv[4][4]={{0.0}};
   PetscReal taua=0.0;
   PetscReal taut=0.0;
   PetscReal taud=0.0;
   PetscReal gij2=0.0;

   taut = 4.0/(dt*dt);


   PetscInt aa,bb,ll;
   for (aa=0;aa<dof;aa++)
           TauTemp[aa][aa] = taut;



   for (aa=0;aa<dof;aa++)
         for (bb=0;bb<dof;bb++)
        	 for (ll=0;ll<dof;ll++){
        		 TauTemp[aa][bb] += G[0][0]*A1_cons[aa][ll]*A1_cons[ll][bb]
                           +  G[0][1]*A1_cons[aa][ll]*A2_cons[ll][bb]
                           +  G[1][0]*A2_cons[aa][ll]*A1_cons[ll][bb]
                           +  G[1][1]*A2_cons[aa][ll]*A2_cons[ll][bb];
         }



    PetscReal m_k=3.0;
    PetscInt k,l;

        for (aa=0;aa<dof;aa++)
               for (bb=0;bb<dof;bb++)
            	   for (ll=0;ll<dof;ll++){


                       TauTemp[aa][bb] += m_k*G[0][0]*K_cons[0][0][aa][ll]*K_cons[0][0][ll][bb]*G[0][0]
                                       +  m_k*G[0][0]*K_cons[0][1][aa][ll]*K_cons[0][0][ll][bb]*G[0][1]
                                       +  m_k*G[0][1]*K_cons[1][0][aa][ll]*K_cons[0][0][ll][bb]*G[0][0]
                                       +  m_k*G[0][1]*K_cons[1][1][aa][ll]*K_cons[0][0][ll][bb]*G[0][1]
                                       +  m_k*G[0][0]*K_cons[0][0][aa][ll]*K_cons[0][1][ll][bb]*G[1][0]
                                       +  m_k*G[0][0]*K_cons[0][1][aa][ll]*K_cons[0][1][ll][bb]*G[1][1]
                                       +  m_k*G[0][1]*K_cons[1][0][aa][ll]*K_cons[0][1][ll][bb]*G[1][0]
                                       +  m_k*G[0][1]*K_cons[1][1][aa][ll]*K_cons[0][1][ll][bb]*G[1][1]

                                       +  m_k*G[1][0]*K_cons[0][0][aa][ll]*K_cons[1][0][ll][bb]*G[0][0]
                                       +  m_k*G[1][0]*K_cons[0][1][aa][ll]*K_cons[1][0][ll][bb]*G[0][1]
                                       +  m_k*G[1][1]*K_cons[1][0][aa][ll]*K_cons[1][0][ll][bb]*G[0][0]
                                       +  m_k*G[1][1]*K_cons[1][1][aa][ll]*K_cons[1][0][ll][bb]*G[0][1]
                                       +  m_k*G[1][0]*K_cons[0][0][aa][ll]*K_cons[1][1][ll][bb]*G[1][0]
                                       +  m_k*G[1][0]*K_cons[0][1][aa][ll]*K_cons[1][1][ll][bb]*G[1][1]
                                       +  m_k*G[1][1]*K_cons[1][0][aa][ll]*K_cons[1][1][ll][bb]*G[1][0]
                                       +  m_k*G[1][1]*K_cons[1][1][aa][ll]*K_cons[1][1][ll][bb]*G[1][1];


               }



  Compute_square_root_inverse(TauTemp, TauTempInv, taut,user);


  //Transform to pressure primitive variables


  for (aa=0;aa<dof;aa++)
        for (bb=0;bb<dof;bb++){
            tau[aa][bb] = A0inv[aa][0]*TauTempInv[0][bb]
                         +A0inv[aa][1]*TauTempInv[1][bb]
                         +A0inv[aa][2]*TauTempInv[2][bb]
                         +A0inv[aa][3]*TauTempInv[3][bb];

       }
 //



  return 0;
}




#undef __FUNCT__
#define __FUNCT__ "ComputeAMatrixConservation"
PetscErrorCode ComputeAMatrixConservation(PetscScalar u[],
        PetscReal (*A0inv)[4],PetscReal (*A1_adv)[4],PetscReal (*A1_p)[4],PetscReal (*A1_pt)[4],PetscReal (*A2_adv)[4],PetscReal (*A2_p)[4],PetscReal (*A2_pt)[4],PetscReal (*K)[2][3][3],PetscReal (*A1_cons)[4],PetscReal (*A2_cons)[4],PetscReal (*K_cons)[2][4][4],void *ctx)
{

    AppCtx *user = (AppCtx *)ctx;

     PetscInt  dof    = user->iga->dof;
     PetscInt i,j,l;

         for (i=0;i<dof;i++){
             for (j=0;j<dof;j++){
                 for (l=0;l<dof;l++){
                        A1_cons[i][j] += (A1_p[i][l]+A1_adv[i][l]+A1_pt[i][l])*A0inv[l][j];
                        A2_cons[i][j] += (A2_p[i][l]+A2_adv[i][l]+A2_pt[i][l])*A0inv[l][j];
                 }
             }
         }

         for (i=0;i<dof-1;i++){
             for (j=0;j<dof;j++){
                 for (l=0;l<dof-1;l++){
                        K_cons[0][0][i+1][j] += K[0][0][i][l]*A0inv[l+1][j];
                        K_cons[0][1][i+1][j] += K[0][1][i][l]*A0inv[l+1][j];
                        K_cons[1][0][i+1][j] += K[1][0][i][l]*A0inv[l+1][j];
                        K_cons[1][1][i+1][j] += K[1][1][i][l]*A0inv[l+1][j];
                 }
             }
         }
  return 0;
}




#undef  __FUNCT__
#define __FUNCT__ "Residual"
PetscErrorCode Residual(IGAPoint pnt,PetscReal dt,
                        PetscReal shift,const PetscScalar *V,
                        PetscReal t,const PetscScalar *U,
                        PetscScalar *Re,void *ctx)
{
	 PetscErrorCode ierr;
   AppCtx *user = (AppCtx *)ctx;

   PetscInt i,j,l, m;
   PetscInt  dof    = pnt->dof;
   PetscInt  dim    = pnt->dim;
   PetscReal lamda  = user->lamda;
   PetscReal mu     = user->mu;
   PetscReal kappa  = user->kappa;
   PetscReal Cv     = user->Cv;

   PetscReal dens0  = 1000.0;
   PetscReal P0     = 100000.0;
   PetscReal B      = 3.31e8;
   PetscReal N      = 7.15;
   PetscReal rhoCR  = dens0*pow((1/B)*(22.02726-P0)+1, 1/N);
   PetscReal Pcr    = 22.02726;

   	 PetscScalar u[dof], u_t[dof];
     IGAPointFormValue(pnt,U,&u[0]);
     IGAPointFormValue(pnt,V,&u_t[0]);

     PetscScalar dens= u[0];
     PetscScalar ux  = u[1];
     PetscScalar uy  = u[2];
     PetscScalar temp= u[3];

     PetscReal P        = Pcr;
     PetscReal fprime   = 0.0;
     PetscReal cs = 0.0;
     if(dens>rhoCR){
     P        = P0+B*(pow(dens/dens0,N))-B;
     fprime   = (1.0/dens)*B*N*(pow(dens/dens0,N));
     cs = (1.0/dens)*B*N*(pow(dens/dens0,N));}

  PetscReal   chi  = lamda+2*mu;
  PetscScalar grad_u[dof][dim];
  //PetscScalar hess_u[dof][dim][dim];
  IGAPointFormGrad(pnt,U,&grad_u[0][0]);
  //IGAPointFormHess(pnt,U,&hess_u[0][0][0]);
  PetscReal InvGradMap[dim][dim];
  IGAPointFormGradMap(pnt,NULL,&InvGradMap[0][0]);

  PetscReal        umi[2]        = {0.0};
  PetscScalar      tau[4][4]     = {{0.0}};
  PetscReal        A1_adv[4][4]  = {{0.0}};
  PetscReal        A1_pt[4][4]   = {{0.0}};
  PetscReal        A2_adv[4][4]  = {{0.0}};
  PetscReal        A2_pt[4][4]   = {{0.0}};
  PetscReal        A1_p[4][4]    = {{0.0}};
  PetscReal        A2_p[4][4]    = {{0.0}};
  PetscReal        A0[4][4]      = {{0.0}};
  PetscReal        A0inv[4][4]   = {{0.0}};
  PetscReal        K[2][2][3][3] = {{{{0.0}}}};
  PetscReal        G[2][2]       = {{0}};

    PetscReal t11,t12,t21,t22;
    t11 = mu*(grad_u[1][0]+grad_u[1][0])+lamda*(grad_u[1][0]+grad_u[2][1]);
    t22 = mu*(grad_u[2][1]+grad_u[2][1])+lamda*(grad_u[1][0]+grad_u[2][1]);
    t12 = mu*(grad_u[1][1]+grad_u[2][0]);
    t21 = mu*(grad_u[2][0]+grad_u[1][1]);

    A0[0][0] = 1.0;
    A0[1][0] = ux;
    A0[1][1] = dens;
    A0[2][0] = uy;
    A0[2][2] = dens;
    A0[3][0] = Cv*temp;
    A0[3][3] = dens*Cv;

    A0inv[0][0] =  1.0;
    A0inv[1][0] = -ux/dens;
    A0inv[1][1] =  1.0/dens;
    A0inv[2][0] = -uy/dens;
    A0inv[2][2] =  1.0/dens;
    A0inv[3][0] =  -temp/dens;
    A0inv[3][3] =  1.0/(Cv*dens);

    //A1
    A1_adv[0][0] = ux;
    A1_adv[0][1] = dens;
    A1_adv[1][0] = ux*ux;
    A1_adv[1][1] = 2.0*dens*ux;
    A1_adv[2][0] = ux*uy;
    A1_adv[2][1] = dens*uy;
    A1_adv[2][2] = dens*ux;
    A1_adv[3][0] = ux*Cv*temp;
    A1_adv[3][1] = dens*Cv*temp;
    A1_adv[3][3] = dens*ux*Cv;

    A1_p[1][0] = fprime;

    A1_pt[3][1] = P-t11;
    A1_pt[3][2] = -t12;

  //A2
  A2_adv[0][0] = uy;
  A2_adv[0][2] = dens;
  A2_adv[1][0] = ux*uy;
  A2_adv[1][1] = dens*uy;
  A2_adv[1][2] = dens*ux;
  A2_adv[2][0] = uy*uy;
  A2_adv[2][2] = 2.0*dens*uy;
  A2_adv[3][0] = uy*Cv*temp;
  A2_adv[3][2] = dens*Cv*temp;
  A2_adv[3][3] = dens*uy*Cv;

  A2_p[2][0] = fprime;

  A2_pt[3][1]=-t21;
  A2_pt[3][2]=P-t22;

  // Viscous terms
  K[0][0][0][0]=chi;
  K[0][0][1][1]=mu;
  K[0][0][2][2]=kappa;

  K[1][1][0][0]=mu;
  K[1][1][1][1]=chi;
  K[1][1][2][2]=kappa;

  K[0][1][1][0]=mu;
  K[0][1][0][1]=lamda;

  K[1][0][0][1]=mu;
  K[1][0][1][0]=lamda;

PetscReal F1[4]={0.0};
PetscReal F2[4]={0.0};
F1[1] = P;//dens*fprime
F2[2] = P;//dens*fprime

for (i=0;i<dim;i++)
for (j=0;j<dim;j++)
for (l=0;l<dim;l++){
G[i][j] += InvGradMap[i][l]*InvGradMap[j][l];
}




PetscReal A1_cons[4][4] = {{0.0}};
PetscReal A2_cons[4][4] = {{0.0}};
PetscReal K_cons[2][2][4][4] = {{{{0.0}}}};

ComputeAMatrixConservation(u,A0inv,A1_adv,A1_p,A1_pt,A2_adv,A2_p,A2_pt,K,A1_cons,A2_cons,K_cons,user);
ierr = DirectTau(G,dt,u,tau,user,A0inv,A1_cons,A2_cons,K_cons,umi);CHKERRQ(ierr);

  PetscReal  *N0 = pnt->shape[0];
  PetscReal (*N1)[dim] = (PetscReal (*)[dim]) pnt->shape[1];

  PetscReal Res[4] = {0.0};

  Res[0]  = -user->F[0];
  Res[1]  = -user->F[1];
  Res[2]  = -user->F[2];
  Res[3]  = -user->F[3];

  for (i=0;i<4;i++){
  for (j=0;j<4;j++){
  	    Res[i] += A0[i][j]*u_t[j] + (A1_adv[i][j]+A1_p[i][j]+A1_pt[i][j])*grad_u[j][0]
         + (A2_adv[i][j]+A2_p[i][j]+A2_pt[i][j])*grad_u[j][1];
  }
  }


  //Stabilization terms DC

   PetscReal hu[3] = {0.0};
   PetscReal A0gradY[4][2]={{0.0}};
   PetscReal tau_m,tau_c,tau_t;

   for (j=0;j<dim;j++)
   for (i=0;i<dof;i++)
   for (l=0;l<dof;l++)
   {
   	A0gradY[i][j] += A0[i][l]*grad_u[l][j];		//ok
   }

   for (i=0;i<dim;i++){
   for (j=0;j<dim;j++){
   	hu[0] += A0gradY[0][i]*G[i][j]*A0gradY[0][j];	//ok
   }
   }

   for (i=0;i<dim;i++){
   for (j=0;j<dim;j++){
   	hu[1] += A0gradY[1][i]*G[i][j]*A0gradY[1][j];   	//ok
   	hu[1] += A0gradY[2][i]*G[i][j]*A0gradY[2][j];		//ok
   }
   }

   for (i=0;i<dim;i++)
   for (j=0;j<dim;j++){
   	hu[2] += A0gradY[3][i]*G[i][j]*A0gradY[3][j];	//ok
   }

   PetscReal Ginv[2][2]={{0.0}};
   PetscReal   det,detinv;

   det  = G[0][0]*G[1][1]-G[0][1]*G[1][0];
   detinv = 1.0/det;
   Ginv[0][0] =  G[1][1]*detinv;
   Ginv[1][0] = -G[1][0]*detinv;
   Ginv[0][1] = -G[0][1]*detinv;
   Ginv[1][1] =  G[0][0]*detinv;

   PetscReal uGu = 0.0;
   for(i=0;i<dim;i++){
     for(j=0;j<dim;j++){
       uGu += u[i+1]*G[i][j]*u[j+1];
     }
   }
   PetscReal k_cap = sqrt(uGu + cs*(G[0][0]+G[1][1]));
   tau_c = (cs*sqrt(Res[0]*Res[0]) + sqrt(ux*ux+uy*uy)*sqrt(Res[1]*Res[1]+Res[2]*Res[2]) + sqrt(Res[3]*Res[3]))/
   (cs*sqrt(hu[0]) + sqrt(ux*ux+uy*uy)*sqrt(hu[1]) + sqrt(hu[2])+1e-6);

   PetscReal DC = 2.0;

   if(DC*tau_c > k_cap){tau_c = k_cap;}else{tau_c = DC*tau_c;}
   tau_m = tau_c;
   tau_t = tau_c;

  PetscScalar (*R)[dof] = (PetscScalar (*)[dof])Re;
  PetscInt a,nen=pnt->nen;
  for (a=0; a<nen; a++) {
    PetscReal Na    = N0[a];
    /* ----- */
    R[a][0]  = -Na*user->F[0];
    R[a][1]  = -Na*user->F[1];
    R[a][2]  = -Na*user->F[2];
    R[a][3]  = -Na*user->F[3];

    for (i=0;i<4;i++){
    for (j=0;j<4;j++){
    	    R[a][i] += (A0[i][j]*u_t[j] + (A1_adv[i][j]+A1_pt[i][j])*grad_u[j][0] + (A2_adv[i][j]+A2_pt[i][j])*grad_u[j][1])*Na;
    }
    }

    for (i=0;i<4;i++){
    	    R[a][i] += -N1[a][0]*F1[i] - N1[a][1]*F2[i];
    }

    for (l=0;l<dim;l++){
    for (m=0;m<dim;m++){
    	for (i=0;i<dof-1;i++){
        for (j=0;j<dof-1;j++){
            		R[a][i+1] += N1[a][l]*K[l][m][i][j]*grad_u[j+1][m];		//ok
        }
    	}
    }
    }

    for (i=0;i<dof;i++)
    for (j=0;j<dof;j++)
    for (l=0;l<dof;l++)
    {
        R[a][i] += ((A1_adv[i][j]+A1_p[i][j]+A1_pt[i][j])*N1[a][0]
         + (A2_adv[i][j]+A2_p[i][j]+A2_pt[i][j])*N1[a][1])*tau[j][l]*Res[l];
    }

    for (i=0;i<dim;i++)
    {
    	    		  R[a][0] += N1[a][i]*tau_c*A0gradY[0][i];
    	    		  R[a][1] += N1[a][i]*tau_m*A0gradY[1][i];
    	    		  R[a][2] += N1[a][i]*tau_m*A0gradY[2][i];
    	    		  R[a][3] += N1[a][i]*tau_t*A0gradY[3][i];
    }
    }
  return 0;
}

#undef  __FUNCT__
#define __FUNCT__ "ResidualFS"
PetscErrorCode ResidualFS(IGAPoint pnt,PetscReal dt,
                        PetscReal shift,const PetscScalar *V,
                        PetscReal t,const PetscScalar *U,
                        PetscScalar *Re,void *ctx)
{
	 PetscErrorCode ierr;
   AppCtx *user = (AppCtx *)ctx;

   PetscInt i,j,l, m;
   PetscInt  dof    = pnt->dof;
   PetscInt  dim    = pnt->dim;
   PetscReal lamda  = user->lamda;
   PetscReal mu     = user->mu;
   PetscReal kappa  = user->kappa;
   PetscReal Cv     = user->Cv;

   PetscReal dens0  = 1000.0;
   PetscReal P0     = 100000.0;
   PetscReal B      = 3.31e8;
   PetscReal N      = 7.15;
   PetscReal rhoCR  = dens0*pow((1/B)*(22.02726-P0)+1, 1/N);
   PetscReal Pcr    = 22.02726;

   	 PetscScalar u[dof], u_t[dof];
     IGAPointFormValue(pnt,U,&u[0]);
     IGAPointFormValue(pnt,V,&u_t[0]);

     PetscScalar dens= u[0];
     PetscScalar ux  = u[1];
     PetscScalar uy  = u[2];
     PetscScalar temp= u[3];

     PetscReal P        = Pcr;
     PetscReal fprime   = 0.0;
     PetscReal cs = 0.0;
     if(dens>rhoCR){
     P        = P0+B*(pow(dens/dens0,N))-B;
     fprime   = (1.0/dens)*B*N*(pow(dens/dens0,N));
     cs = (1.0/dens)*B*N*(pow(dens/dens0,N));}

  PetscReal   chi  = lamda+2*mu;
  PetscScalar grad_u[dof][dim];
  //PetscScalar hess_u[dof][dim][dim];
  IGAPointFormGrad (pnt,U,&grad_u[0][0]);
  //IGAPointFormHess(pnt,U,&hess_u[0][0][0]);
  PetscReal InvGradMap[dim][dim];
  IGAPointFormGradMap(pnt,NULL,&InvGradMap[0][0]);

  PetscReal        umi[2]        = {0.0};
  PetscScalar      tau[4][4]     = {{0.0}};
  PetscReal        A1_adv[4][4]  = {{0.0}};
  PetscReal        A1_pt[4][4]   = {{0.0}};
  PetscReal        A2_adv[4][4]  = {{0.0}};
  PetscReal        A2_pt[4][4]   = {{0.0}};
  PetscReal        A1_p[4][4]    = {{0.0}};
  PetscReal        A2_p[4][4]    = {{0.0}};
  PetscReal        A0[4][4]      = {{0.0}};
  PetscReal        A0inv[4][4]   = {{0.0}};
  PetscReal        K[2][2][3][3] = {{{{0.0}}}};
  PetscReal        G[2][2]       = {{0}};

    PetscReal t11,t12,t21,t22;
    t11 = mu*(grad_u[1][0]+grad_u[1][0])+lamda*(grad_u[1][0]+grad_u[2][1]);
    t22 = mu*(grad_u[2][1]+grad_u[2][1])+lamda*(grad_u[1][0]+grad_u[2][1]);
    t12 = mu*(grad_u[1][1]+grad_u[2][0]);
    t21 = mu*(grad_u[2][0]+grad_u[1][1]);

    A0[0][0] = 1.0;
    A0[1][0] = ux;
    A0[1][1] = dens;
    A0[2][0] = uy;
    A0[2][2] = dens;
    A0[3][0] = Cv*temp;
    A0[3][3] = dens*Cv;

    A0inv[0][0] =  1.0;
    A0inv[1][0] = -ux/dens;
    A0inv[1][1] =  1.0/dens;
    A0inv[2][0] = -uy/dens;
    A0inv[2][2] =  1.0/dens;
    A0inv[3][0] =  -temp/dens;
    A0inv[3][3] =  1.0/(Cv*dens);

    //A1
    A1_adv[0][0] = ux;
    A1_adv[0][1] = dens;
    A1_adv[1][0] = ux*ux;
    A1_adv[1][1] = 2.0*dens*ux;
    A1_adv[2][0] = ux*uy;
    A1_adv[2][1] = dens*uy;
    A1_adv[2][2] = dens*ux;
    A1_adv[3][0] = ux*Cv*temp;
    A1_adv[3][1] = dens*Cv*temp;
    A1_adv[3][3] = dens*ux*Cv;

    A1_p[1][0] = fprime;

    A1_pt[3][1] = P-t11;
    A1_pt[3][2] = -t12;

  //A2
  A2_adv[0][0] = uy;
  A2_adv[0][2] = dens;
  A2_adv[1][0] = ux*uy;
  A2_adv[1][1] = dens*uy;
  A2_adv[1][2] = dens*ux;
  A2_adv[2][0] = uy*uy;
  A2_adv[2][2] = 2.0*dens*uy;
  A2_adv[3][0] = uy*Cv*temp;
  A2_adv[3][2] = dens*Cv*temp;
  A2_adv[3][3] = dens*uy*Cv;

  A2_p[2][0] = fprime;

  A2_pt[3][1]=-t21;
  A2_pt[3][2]=P-t22;

  // Viscous terms
  K[0][0][0][0]=chi;
  K[0][0][1][1]=mu;
  K[0][0][2][2]=kappa;

  K[1][1][0][0]=mu;
  K[1][1][1][1]=chi;
  K[1][1][2][2]=kappa;

  K[0][1][1][0]=mu;
  K[0][1][0][1]=lamda;

  K[1][0][0][1]=mu;
  K[1][0][1][0]=lamda;

PetscReal F1[4]={0.0};
PetscReal F2[4]={0.0};
F1[1] = P;//dens*fprime
F2[2] = P;//dens*fprime

for (i=0;i<dim;i++)
for (j=0;j<dim;j++)
for (l=0;l<dim;l++){
G[i][j] += InvGradMap[i][l]*InvGradMap[j][l];
}

PetscReal A1_cons[4][4] = {{0.0}};
PetscReal A2_cons[4][4] = {{0.0}};
PetscReal K_cons[2][2][4][4] = {{{{0.0}}}};

ComputeAMatrixConservation(u,A0inv,A1_adv,A1_p,A1_pt,A2_adv,A2_p,A2_pt,K,A1_cons,A2_cons,K_cons,user);
ierr = DirectTau(G,dt,u,tau,user,A0inv,A1_cons,A2_cons,K_cons,umi);CHKERRQ(ierr);

  PetscReal  *N0 = pnt->shape[0];
  PetscReal (*N1)[dim] = (PetscReal (*)[dim]) pnt->shape[1];

  PetscReal Res[4] = {0.0};

  Res[0]  = -user->F[0];
  Res[1]  = -user->F[1];
  Res[2]  = -user->F[2];
  Res[3]  = -user->F[3];

  for (i=0;i<4;i++){
  for (j=0;j<4;j++){
  	    Res[i] += A0[i][j]*u_t[j] + (A1_adv[i][j]+A1_p[i][j]+A1_pt[i][j])*grad_u[j][0]
         + (A2_adv[i][j]+A2_p[i][j]+A2_pt[i][j])*grad_u[j][1];
  }
  }


  //Stabilization terms DC

   PetscReal hu[3] = {0.0};
   PetscReal A0gradY[4][2]={{0.0}};
   PetscReal tau_m,tau_c,tau_t;

   for (j=0;j<dim;j++)
   for (i=0;i<dof;i++)
   for (l=0;l<dof;l++)
   {
   	A0gradY[i][j] += A0[i][l]*grad_u[l][j];		//ok
   }

   for (i=0;i<dim;i++){
   for (j=0;j<dim;j++){
   	hu[0] += A0gradY[0][i]*G[i][j]*A0gradY[0][j];	//ok
   }
   }
   /*if (hu > 1e-12){
   tau_c = sqrt(Res[0]*Res[0]/hu);
   }else{
   tau_c = 0.0;
 }*/


   for (i=0;i<dim;i++){
   for (j=0;j<dim;j++){
   	hu[1] += A0gradY[1][i]*G[i][j]*A0gradY[1][j];   	//ok
   	hu[1] += A0gradY[2][i]*G[i][j]*A0gradY[2][j];		//ok
   }
   }
   /*if (hu > 1e-12){
   tau_m = sqrt( (Res[1]*Res[1] + Res[2]*Res[2]) /hu);
   }else{
   tau_m = 0.0;
 }*/


   for (i=0;i<dim;i++)
   for (j=0;j<dim;j++){
   	hu[2] += A0gradY[3][i]*G[i][j]*A0gradY[3][j];	//ok
   }

   /*if (hu > 1e-12){
   tau_t = sqrt(Res[3]*Res[3]/hu);
   }else{
   tau_t = 0.0;
 }*/

   PetscReal Ginv[2][2]={{0.0}};
   PetscReal   det,detinv;

   det  = G[0][0]*G[1][1]-G[0][1]*G[1][0];
   detinv = 1.0/det;
   Ginv[0][0] =  G[1][1]*detinv;
   Ginv[1][0] = -G[1][0]*detinv;
   Ginv[0][1] = -G[0][1]*detinv;
   Ginv[1][1] =  G[0][0]*detinv;

   PetscReal uGu = 0.0;
   for(i=0;i<dim;i++){
     for(j=0;j<dim;j++){
       uGu += u[i+1]*G[i][j]*u[j+1];
     }
   }
   PetscReal k_cap = sqrt(uGu + cs*(G[0][0]+G[1][1]));
   tau_c = (cs*sqrt(Res[0]*Res[0]) + sqrt(ux*ux+uy*uy)*sqrt(Res[1]*Res[1]+Res[2]*Res[2]) + sqrt(Res[3]*Res[3]))/
   (cs*sqrt(hu[0]) + sqrt(ux*ux+uy*uy)*sqrt(hu[1]) + sqrt(hu[2])+1e-6);

   PetscReal DC = 2.0;

   if(DC*tau_c > k_cap){tau_c = k_cap;}else{tau_c = DC*tau_c;}
   tau_m = tau_c;
   tau_t = tau_c;

  PetscScalar (*R)[dof] = (PetscScalar (*)[dof])Re;
  PetscInt a,nen=pnt->nen;
  for (a=0; a<nen; a++) {
    PetscReal Na    = N0[a];
    /* ----- */
    R[a][0]  = -Na*user->F[0];
    R[a][1]  = -Na*user->F[1];
    R[a][2]  = -Na*user->F[2];
    R[a][3]  = -Na*user->F[3];

    for (i=0;i<4;i++){
    for (j=0;j<4;j++){
    	    R[a][i] += (A0[i][j]*u_t[j] + (A1_adv[i][j]+A1_pt[i][j])*grad_u[j][0] + (A2_adv[i][j]+A2_pt[i][j])*grad_u[j][1])*Na;
    	    //R[a][i] += -N1[a][0]*A1_p[i][j]*u[j] - N1[a][1]*A2_p[i][j]*u[j];
		      //R[a][i] += (A0[i][j]*u_t[j] + A1_c[i][j]*grad_u[j][0] + A2_c[i][j]*grad_u[j][1])*Na;
    }
    }

    for (i=0;i<4;i++){
    	    R[a][i] += -N1[a][0]*F1[i] - N1[a][1]*F2[i];
    }


    for (l=0;l<dim;l++){
    for (m=0;m<dim;m++){
    	for (i=0;i<dof-1;i++){
        for (j=0;j<dof-1;j++){
            		R[a][i+1] += N1[a][l]*K[l][m][i][j]*grad_u[j+1][m];		//ok
        }
    	}
    }
    }

    for (i=0;i<dof;i++)
    for (j=0;j<dof;j++)
    for (l=0;l<dof;l++)
    {
        R[a][i] += ((A1_adv[i][j]+A1_p[i][j]+A1_pt[i][j])*N1[a][0]
         + (A2_adv[i][j]+A2_p[i][j]+A2_pt[i][j])*N1[a][1])*tau[j][l]*Res[l];
    }

    for (i=0;i<dim;i++)
    {
    	    		  R[a][0] += N1[a][i]*tau_c*A0gradY[0][i];
    	    		  R[a][1] += N1[a][i]*tau_m*A0gradY[1][i];
    	    		  R[a][2] += N1[a][i]*tau_m*A0gradY[2][i];
    	    		  R[a][3] += N1[a][i]*tau_t*A0gradY[3][i];
    }
    }
  return 0;
}

#undef  __FUNCT__
#define __FUNCT__ "ResidualRDX"
PetscErrorCode ResidualRDX(IGAPoint pnt,PetscReal dt,
                        PetscReal shift,const PetscScalar *V,
                        PetscReal t,const PetscScalar *U,
                        PetscScalar *Re,void *ctx)
{
	 PetscErrorCode ierr;
   AppCtx *user = (AppCtx *)ctx;

   PetscInt i,j,l, m;
   PetscInt  dof    = pnt->dof;
   PetscInt  dim    = pnt->dim;
   PetscReal lamda  = user->lamda;
   PetscReal mu     = user->mu;
   PetscReal kappa  = user->kappa;
   PetscReal Cv     = user->Cv;

   PetscReal dens0  = 1659.0;
   PetscReal P0     = 100000.0;
   PetscReal A      = 495.1e9;
   PetscReal B      = 7.21e9;
   PetscReal C      = 1.62e9;
   PetscReal R1     = 4.387;
   PetscReal R2     = 0.9954;
   PetscReal omega  = 0.3469;
   PetscReal E0     = 5877.9e3;
   PetscReal Pcr    = 2.0e11;

   	 PetscScalar u[dof], u_t[dof];
     IGAPointFormValue(pnt,U,&u[0]);
     IGAPointFormValue(pnt,V,&u_t[0]);
     PetscScalar dens= u[0];
     PetscScalar ux  = u[1];
     PetscScalar uy  = u[2];
     PetscScalar temp= u[3];

     PetscReal nu     = dens/dens0;
     PetscReal Ptest  = A*exp(-R1*nu) + B*exp(-R2*nu) + C/(pow(nu, 1+omega)); //We can consider just one branch of the JWL
     //Equation, as they intersect at 2E11, to determine if P>2E11. Then, we will choose the right function based
     //on this value.
     PetscReal P        = Pcr;
     PetscReal fprime   = 0.0;
     PetscReal cs       = 0.0;

     if(Ptest > Pcr){
       P        = Ptest;
       fprime   = A*exp(-R1*nu)*R1*nu/dens +  B*exp(-R2*nu)*R2*nu/dens  +  C*(1+omega)*pow(dens,omega)/pow(dens0,1+omega);
       cs       = fprime;
     }else{
       P        = A*(1-omega/(R1*nu))*exp(-R1*nu) + B*(1-omega/(R2*nu))*exp(-R2*nu) + omega*dens0*E0/nu;
       fprime   = ((A*omega/(R1*dens0))*(R1*nu-1)+nu*A*R1/dens)*exp(-R1*nu) + ((B*omega/(R2*dens0))*(R2*nu-1)+nu*B*R2/dens)*exp(-R2*nu) + omega*E0;
       cs       = fprime;
     }

  PetscReal   chi  = lamda+2*mu;
  PetscScalar grad_u[dof][dim];
  //PetscScalar hess_u[dof][dim][dim];
  IGAPointFormGrad (pnt,U,&grad_u[0][0]);
  //IGAPointFormHess(pnt,U,&hess_u[0][0][0]);
  PetscReal InvGradMap[dim][dim];
  IGAPointFormGradMap(pnt,NULL,&InvGradMap[0][0]);

  PetscReal        umi[2]        = {0.0};
  PetscScalar      tau[4][4]     = {{0.0}};
  PetscReal        A1_adv[4][4]  = {{0.0}};
  PetscReal        A1_pt[4][4]   = {{0.0}};
  PetscReal        A2_adv[4][4]  = {{0.0}};
  PetscReal        A2_pt[4][4]   = {{0.0}};
  PetscReal        A1_p[4][4]    = {{0.0}};
  PetscReal        A2_p[4][4]    = {{0.0}};
  PetscReal        A0[4][4]      = {{0.0}};
  PetscReal        A0inv[4][4]   = {{0.0}};
  PetscReal        K[2][2][3][3] = {{{{0.0}}}};
  PetscReal        G[2][2]       = {{0}};

    PetscReal t11,t12,t21,t22;
    t11 = mu*(grad_u[1][0]+grad_u[1][0])+lamda*(grad_u[1][0]+grad_u[2][1]);
    t22 = mu*(grad_u[2][1]+grad_u[2][1])+lamda*(grad_u[1][0]+grad_u[2][1]);
    t12 = mu*(grad_u[1][1]+grad_u[2][0]);
    t21 = mu*(grad_u[2][0]+grad_u[1][1]);

    A0[0][0] = 1.0;
    A0[1][0] = ux;
    A0[1][1] = dens;
    A0[2][0] = uy;
    A0[2][2] = dens;
    A0[3][0] = Cv*temp;
    A0[3][3] = dens*Cv;

    A0inv[0][0] =  1.0;
    A0inv[1][0] = -ux/dens;
    A0inv[1][1] =  1.0/dens;
    A0inv[2][0] = -uy/dens;
    A0inv[2][2] =  1.0/dens;
    A0inv[3][0] =  -temp/dens;
    A0inv[3][3] =  1.0/(Cv*dens);

    //A1
    A1_adv[0][0] = ux;
    A1_adv[0][1] = dens;
    A1_adv[1][0] = ux*ux;
    A1_adv[1][1] = 2.0*dens*ux;
    A1_adv[2][0] = ux*uy;
    A1_adv[2][1] = dens*uy;
    A1_adv[2][2] = dens*ux;
    A1_adv[3][0] = ux*Cv*temp;
    A1_adv[3][1] = dens*Cv*temp;
    A1_adv[3][3] = dens*ux*Cv;

    A1_p[1][0] = fprime;

    A1_pt[3][1] = P-t11;
    A1_pt[3][2] = -t12;

  //A2
  A2_adv[0][0] = uy;
  A2_adv[0][2] = dens;
  A2_adv[1][0] = ux*uy;
  A2_adv[1][1] = dens*uy;
  A2_adv[1][2] = dens*ux;
  A2_adv[2][0] = uy*uy;
  A2_adv[2][2] = 2.0*dens*uy;
  A2_adv[3][0] = uy*Cv*temp;
  A2_adv[3][2] = dens*Cv*temp;
  A2_adv[3][3] = dens*uy*Cv;

  A2_p[2][0] = fprime;

  A2_pt[3][1]=-t21;
  A2_pt[3][2]=P-t22;

  // Viscous terms
  K[0][0][0][0]=chi;
  K[0][0][1][1]=mu;
  K[0][0][2][2]=kappa;

  K[1][1][0][0]=mu;
  K[1][1][1][1]=chi;
  K[1][1][2][2]=kappa;

  K[0][1][1][0]=mu;
  K[0][1][0][1]=lamda;

  K[1][0][0][1]=mu;
  K[1][0][1][0]=lamda;

PetscReal F1[4]={0.0};
PetscReal F2[4]={0.0};
F1[1] = P;//dens*fprime
F2[2] = P;//dens*fprime

for (i=0;i<dim;i++)
for (j=0;j<dim;j++)
for (l=0;l<dim;l++){
G[i][j] += InvGradMap[i][l]*InvGradMap[j][l];
}

PetscReal A1_cons[4][4] = {{0.0}};
PetscReal A2_cons[4][4] = {{0.0}};
PetscReal K_cons[2][2][4][4] = {{{{0.0}}}};

ComputeAMatrixConservation(u,A0inv,A1_adv,A1_p,A1_pt,A2_adv,A2_p,A2_pt,K,A1_cons,A2_cons,K_cons,user);
ierr = DirectTau(G,dt,u,tau,user,A0inv,A1_cons,A2_cons,K_cons,umi);CHKERRQ(ierr);

  PetscReal  *N0 = pnt->shape[0];
  PetscReal (*N1)[dim] = (PetscReal (*)[dim]) pnt->shape[1];

  PetscReal Res[4] = {0.0};

  Res[0]  = -user->F[0];
  Res[1]  = -user->F[1];
  Res[2]  = -user->F[2];
  Res[3]  = -user->F[3];

  for (i=0;i<4;i++){
  for (j=0;j<4;j++){
  	    Res[i] += A0[i][j]*u_t[j] + (A1_adv[i][j]+A1_p[i][j]+A1_pt[i][j])*grad_u[j][0]
         + (A2_adv[i][j]+A2_p[i][j]+A2_pt[i][j])*grad_u[j][1];
  }
  }


  //Stabilization terms DC

   PetscReal hu[3] = {0.0};
   PetscReal A0gradY[4][2]={{0.0}};
   PetscReal tau_m,tau_c,tau_t;

   for (j=0;j<dim;j++)
   for (i=0;i<dof;i++)
   for (l=0;l<dof;l++)
   {
   	A0gradY[i][j] += A0[i][l]*grad_u[l][j];		//ok
   }

   for (i=0;i<dim;i++){
   for (j=0;j<dim;j++){
   	hu[0] += A0gradY[0][i]*G[i][j]*A0gradY[0][j];	//ok
   }
   }
   /*if (hu > 1e-12){
   tau_c = sqrt(Res[0]*Res[0]/hu);
   }else{
   tau_c = 0.0;
 }*/


   for (i=0;i<dim;i++){
   for (j=0;j<dim;j++){
   	hu[1] += A0gradY[1][i]*G[i][j]*A0gradY[1][j];   	//ok
   	hu[1] += A0gradY[2][i]*G[i][j]*A0gradY[2][j];		//ok
   }
   }
   /*if (hu > 1e-12){
   tau_m = sqrt( (Res[1]*Res[1] + Res[2]*Res[2]) /hu);
   }else{
   tau_m = 0.0;
 }*/


   for (i=0;i<dim;i++)
   for (j=0;j<dim;j++){
   	hu[2] += A0gradY[3][i]*G[i][j]*A0gradY[3][j];	//ok
   }

   /*if (hu > 1e-12){
   tau_t = sqrt(Res[3]*Res[3]/hu);
   }else{
   tau_t = 0.0;
 }*/

   PetscReal Ginv[2][2]={{0.0}};
   PetscReal   det,detinv;

   det  = G[0][0]*G[1][1]-G[0][1]*G[1][0];
   detinv = 1.0/det;
   Ginv[0][0] =  G[1][1]*detinv;
   Ginv[1][0] = -G[1][0]*detinv;
   Ginv[0][1] = -G[0][1]*detinv;
   Ginv[1][1] =  G[0][0]*detinv;

   PetscReal uGu = 0.0;
   for(i=0;i<dim;i++){
     for(j=0;j<dim;j++){
       uGu += u[i+1]*G[i][j]*u[j+1];
     }
   }
   PetscReal k_cap = sqrt(uGu + cs*(G[0][0]+G[1][1]));
   tau_c = (cs*sqrt(Res[0]*Res[0]) + sqrt(ux*ux+uy*uy)*sqrt(Res[1]*Res[1]+Res[2]*Res[2]) + sqrt(Res[3]*Res[3]))/
   (cs*sqrt(hu[0]) + sqrt(ux*ux+uy*uy)*sqrt(hu[1]) + sqrt(hu[2])+1e-6);

   PetscReal DC = 2.0;

   if(DC*tau_c > k_cap){tau_c = k_cap;}else{tau_c = DC*tau_c;}
   tau_m = tau_c;
   tau_t = tau_c;

  PetscScalar (*R)[dof] = (PetscScalar (*)[dof])Re;
  PetscInt a,nen=pnt->nen;
  for (a=0; a<nen; a++) {
    PetscReal Na    = N0[a];
    /* ----- */
    R[a][0]  = -Na*user->F[0];
    R[a][1]  = -Na*user->F[1];
    R[a][2]  = -Na*user->F[2];
    R[a][3]  = -Na*user->F[3];

    for (i=0;i<4;i++){
    for (j=0;j<4;j++){
    	    R[a][i] += (A0[i][j]*u_t[j] + (A1_adv[i][j]+A1_pt[i][j])*grad_u[j][0] + (A2_adv[i][j]+A2_pt[i][j])*grad_u[j][1])*Na;
    	    //R[a][i] += -N1[a][0]*A1_p[i][j]*u[j] - N1[a][1]*A2_p[i][j]*u[j];
		      //R[a][i] += (A0[i][j]*u_t[j] + A1_c[i][j]*grad_u[j][0] + A2_c[i][j]*grad_u[j][1])*Na;
    }
    }

    for (i=0;i<4;i++){
    	    R[a][i] += -N1[a][0]*F1[i] - N1[a][1]*F2[i];
    }


    for (l=0;l<dim;l++){
    for (m=0;m<dim;m++){
    	for (i=0;i<dof-1;i++){
        for (j=0;j<dof-1;j++){
            		R[a][i+1] += N1[a][l]*K[l][m][i][j]*grad_u[j+1][m];		//ok
        }
    	}
    }
    }

    for (i=0;i<dof;i++)
    for (j=0;j<dof;j++)
    for (l=0;l<dof;l++)
    {
        R[a][i] += ((A1_adv[i][j]+A1_p[i][j]+A1_pt[i][j])*N1[a][0]
         + (A2_adv[i][j]+A2_p[i][j]+A2_pt[i][j])*N1[a][1])*tau[j][l]*Res[l];
    }

    for (i=0;i<dim;i++)
    {
    	    		  R[a][0] += N1[a][i]*tau_c*A0gradY[0][i];
    	    		  R[a][1] += N1[a][i]*tau_m*A0gradY[1][i];
    	    		  R[a][2] += N1[a][i]*tau_m*A0gradY[2][i];
    	    		  R[a][3] += N1[a][i]*tau_t*A0gradY[3][i];
    }
    }
  return 0;
}

#undef  __FUNCT__
#define __FUNCT__ "ResidualBound"
PetscErrorCode ResidualBound(IGAPoint pnt,PetscReal dt,
                        PetscReal shift,const PetscScalar *V,
                        PetscReal t,const PetscScalar *U,
                        PetscScalar *Re,void *ctx)
{

	 PetscErrorCode ierr;
   AppCtx *user = (AppCtx *)ctx;

   PetscInt  i,j,l;
   PetscInt  dof    = pnt->dof;
   PetscInt  dim    = pnt->dim;
   PetscReal lamda  = user->lamda;
   PetscReal mu     = user->mu;
   PetscReal kappa  = user->kappa;
   PetscReal Cv     = user->Cv;
   PetscReal RR     = user->R;

   PetscReal dens0  = 1000.0;
   PetscReal P0     = 100000.0;
   PetscReal B      = 3.31e8;
   PetscReal N      = 7.15;
   PetscReal rhoCR  = dens0*pow((1/B)*(22.02726-P0)+1, 1/N);
   PetscReal Pcr    = 22.02726;

     PetscScalar u[dof];
     IGAPointFormValue(pnt,U,&u[0]);


     PetscScalar dens = u[0];
     PetscScalar ux   = u[1];
     PetscScalar uy   = u[2];
     PetscScalar temp = u[3];

     PetscReal P        = Pcr;
     PetscReal fprime   = 0.0;

     if(dens>rhoCR){
     P        = P0+B*(pow(dens/dens0,N))-B;
     fprime   = (1.0/dens)*B*N*(pow(dens/dens0,N));}
     // Determine if control point is Explosive gas products or if it is
     // water. This will determine the constituative model used to enforce
     // weakly pressure conditions.
    if(sqrt((dens-1000.0)*(dens-1000.0))>0.001){
       PetscReal dens0  = 1659.0;
       PetscReal P0     = 100000.0;
       PetscReal A      = 495.1e9;
       PetscReal B      = 7.21e9;
       PetscReal C      = 1.62e9;
       PetscReal R1     = 4.387;
       PetscReal R2     = 0.9954;
       PetscReal omega  = 0.3469;
       PetscReal E0     = 5877.9e3;
       PetscReal Pcr    = 2.0e11;
       Cv = 143.3;
       PetscReal nu     = dens/dens0;
       PetscReal Ptest  = A*exp(-R1*nu) + B*exp(-R2*nu) + C/(pow(nu, 1+omega)); //We can consider just one branch of the JWL
       //Equation, as they intersect at 2E11, to determine if P>2E11. Then, we will choose the right function based
       //on this value.
       PetscReal P        = Pcr;
       PetscReal fprime   = 0.0;
       PetscReal cs       = 0.0;

       if(Ptest > Pcr){
         P        = Ptest;
         fprime   = A*exp(-R1*nu)*R1*nu/dens +  B*exp(-R2*nu)*R2*nu/dens  +  C*(1+omega)*pow(dens,omega)/pow(dens0,1+omega);
         cs       = fprime;
       }else{
         P        = A*(1-omega/(R1*nu))*exp(-R1*nu) + B*(1-omega/(R2*nu))*exp(-R2*nu) + omega*dens0*E0/nu;
         fprime   = ((A*omega/(R1*dens0))*(R1*nu-1)+nu*A*R1/dens)*exp(-R1*nu) + ((B*omega/(R2*dens0))*(R2*nu-1)+nu*B*R2/dens)*exp(-R2*nu) + omega*E0;
         cs       = fprime;
       }
     }

  PetscReal      A0[4][4] = {{0.0}};

  A0[0][0] = 1.0;
  A0[1][0] = ux;
  A0[1][1] = dens;
  A0[2][0] = uy;
  A0[2][2] = dens;
  A0[3][0] = Cv*temp;
  A0[3][3] = dens*Cv;

  PetscReal F1[4]={0.0};
  PetscReal F2[4]={0.0};
  F1[1] = P;
  F2[2] = P;

  PetscReal  *N0 = pnt->shape[0];
  PetscReal  *norm;
  norm       = pnt->normal;

  PetscScalar (*R)[dof] = (PetscScalar (*)[dof])Re;
  PetscInt a,nen=pnt->nen;
  for (a=0; a<nen; a++) {
    PetscReal Na    = N0[a];
    for (i=0;i<4;i++){
		    R[a][i] += Na*F1[i]*norm[0] + Na*F2[i]*norm[1];
    }
  }
  return 0;
}



#undef  __FUNCT__
#define __FUNCT__ "Tangent"
PetscErrorCode Tangent(IGAPoint pnt,PetscReal dt,
                       PetscReal shift,const PetscScalar *V,
                       PetscReal t,const PetscScalar *U,
                       PetscScalar *Ke,void *ctx)
{
  AppCtx *user = (AppCtx *)ctx;

  PetscInt i,j,l;
  PetscInt  dof    = pnt->dof;
  PetscInt  dim    = pnt->dim;
  PetscReal Cv     = user->Cv;

    PetscScalar u[dof];
    IGAPointFormValue(pnt,U,&u[0]);
    PetscReal *N0 = pnt->shape[0];

    PetscScalar dens= u[0];
    PetscScalar ux  = u[1];
    PetscScalar uy  = u[2];
    PetscScalar temp= u[3];

  PetscReal A0[4][4]={{0.0}};

  A0[0][0] = 1.0;
  A0[1][0] = ux;
  A0[1][1] = dens;
  A0[2][0] = uy;
  A0[2][2] = dens;
  A0[3][0] = Cv*temp;
  A0[3][3] = dens*Cv;

  PetscInt a,nen=pnt->nen;

  PetscScalar (*KK)[dof][nen][dof] = (PetscScalar (*)[dof][nen][dof])Ke;
  for (a=0; a<nen; a++) {
      KK[a][0][a][0] = A0[0][0]*shift*N0[a];
      KK[a][0][a][1] = A0[0][1]*shift*N0[a];
      KK[a][0][a][2] = A0[0][2]*shift*N0[a];
      KK[a][0][a][3] = A0[0][3]*shift*N0[a];
      KK[a][1][a][0] = A0[1][0]*shift*N0[a];
      KK[a][1][a][1] = A0[1][1]*shift*N0[a];
      KK[a][1][a][2] = A0[1][2]*shift*N0[a];
      KK[a][1][a][3] = A0[1][3]*shift*N0[a];
      KK[a][2][a][0] = A0[2][0]*shift*N0[a];
      KK[a][2][a][1] = A0[2][1]*shift*N0[a];
      KK[a][2][a][2] = A0[2][2]*shift*N0[a];
      KK[a][2][a][3] = A0[2][3]*shift*N0[a];
      KK[a][3][a][0] = A0[3][0]*shift*N0[a];
      KK[a][3][a][1] = A0[3][1]*shift*N0[a];
      KK[a][3][a][2] = A0[3][2]*shift*N0[a];
      KK[a][3][a][3] = A0[3][3]*shift*N0[a];
  }
  return 0;
}



#undef  __FUNCT__
#define __FUNCT__ "IGAComputeIJacobianComp"
PetscErrorCode IGAComputeIJacobianComp(IGA iga,PetscReal dt,
                                   PetscReal a,Vec vecV,
                                   PetscReal t,Vec vecU,
                                   Mat matJ,AppCtx *user)
{
  Vec               localV;
  Vec               localU;
  const PetscScalar *arrayV;
  const PetscScalar *arrayU;
  IGAElement        element;
  IGAPoint          point;
  IGAUserIJacobian  IJacobian;
  void              *ctx;
  PetscScalar       *V,*U,*J,*K;
  PetscErrorCode    ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidHeaderSpecific(vecV,VEC_CLASSID,4);
  PetscValidHeaderSpecific(vecU,VEC_CLASSID,6);
  PetscValidHeaderSpecific(matJ,MAT_CLASSID,7);
  IGACheckSetUp(iga,1);
  IGACheckUserOp(iga,1,IJacobian);

  /* Clear global matrix J*/
  ierr = MatZeroEntries(matJ);CHKERRQ(ierr);

  /* Get local vectors V,U and arrays */
  ierr = IGAGetLocalVecArray(iga,vecV,&localV,&arrayV);CHKERRQ(ierr);
  ierr = IGAGetLocalVecArray(iga,vecU,&localU,&arrayU);CHKERRQ(ierr);
  ierr = PetscLogEventBegin(IGA_FormIJacobian,iga,vecV,vecU,matJ);CHKERRQ(ierr);

  /* Element Loop */
  ierr = IGABeginElement(iga,&element);CHKERRQ(ierr);
  while (IGANextElement(iga,element)) {
    ierr = IGAElementGetWorkVal(element,&V);CHKERRQ(ierr);
    ierr = IGAElementGetWorkVal(element,&U);CHKERRQ(ierr);
    ierr = IGAElementGetWorkMat(element,&J);CHKERRQ(ierr);
    ierr = IGAElementGetValues(element,arrayV,V);CHKERRQ(ierr);
    ierr = IGAElementGetValues(element,arrayU,U);CHKERRQ(ierr);
    ierr = IGAElementFixValues(element,U);CHKERRQ(ierr);
    /* UserIJacobian loop */
//    while (IGAElementNextUserIJacobian(element,&IJacobian,&ctx)) {
      /* Quadrature loop */
      ierr = IGAElementBeginPoint(element,&point);CHKERRQ(ierr);
      while (IGAElementNextPoint(element,point)) {
        ierr = IGAPointGetWorkMat(point,&K);CHKERRQ(ierr);
        ierr = Tangent(point,dt,a,V,t,U,K,user);CHKERRQ(ierr);
        ierr = IGAPointAddMat(point,K,J);CHKERRQ(ierr);
      }
      ierr = IGAElementEndPoint(element,&point);CHKERRQ(ierr);
    //}
    ierr = IGAElementFixJacobian(element,J);CHKERRQ(ierr);
    ierr = IGAElementAssembleMat(element,J,matJ);CHKERRQ(ierr);
  }
  ierr = IGAEndElement(iga,&element);CHKERRQ(ierr);

  ierr = PetscLogEventEnd(IGA_FormIJacobian,iga,vecV,vecU,matJ);CHKERRQ(ierr);

  /* Get local vectors V,U and arrays */
  ierr = IGARestoreLocalVecArray(iga,vecV,&localV,&arrayV);CHKERRQ(ierr);
  ierr = IGARestoreLocalVecArray(iga,vecU,&localU,&arrayU);CHKERRQ(ierr);

  /* Assemble global matrix J*/
  ierr = MatAssemblyBegin(matJ,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd  (matJ,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr = IGAElementDestroy(&element);CHKERRQ(ierr);
  ierr = IGAPointDestroy(&point);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "outputTXT"
PetscErrorCode outputTXT (PARAMETERS *par)
{

  PetscInt count=par->numNodes;
  PetscInt outCount=par->stepNumber;
  PetscInt counter=par->stepNumber/par->FreqResults;


  if (outCount % par->FreqResults == 0){


  // ##################################################
  //                  Geometry File
  // ##################################################


  FILE *fileResult;
  char filenameResults[50];

  sprintf(filenameResults,"./Meshless.%d.geo",counter);


  fileResult=fopen(filenameResults,"wt");
  if (fileResult == NULL)
      {
  	printf("Error opening file!\n");
  	exit(1);
      }

  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"Meshless \n");
  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"node \n");
  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"node id given \n");
  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"element id given \n");
  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"coordinates \n");


  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"%d \n", count);

  PetscInt i,j;

  for(i=0; i<count;i++)
  {
    PetscFPrintf(PETSC_COMM_WORLD,fileResult,"%d   %.4e   %.4e   %.4e \n", i+1, par->puntos[i].currentCoord[0], par->puntos[i].currentCoord[1], 0.0);
  }

  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"part     1\n");
  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"todo \n");
  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"point \n");
  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"%d \n",count);

  for(j=1;j<=count;j++){
  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"%d  %d \n",j,j);
  }

  fclose(fileResult);



  // ##################################################
   //                  Case File
   // ##################################################

//  PetscPrintf(PETSC_COMM_WORLD," outCount %d \n",outCount );

  if (outCount == 0){

  sprintf(filenameResults,"./Meshfree.case");
  fileResult=fopen(filenameResults,"wt");

   if (fileResult == NULL)
       {
   	printf("Error opening file!\n");
   	exit(1);
       }
   PetscInt numSteps = (int) (par->finalTime/par->timeStep);

   PetscFPrintf(PETSC_COMM_WORLD,fileResult,"#BOF: meshless.case \n");
   PetscFPrintf(PETSC_COMM_WORLD,fileResult,"\n");
   PetscFPrintf(PETSC_COMM_WORLD,fileResult,"FORMAT \n");
   PetscFPrintf(PETSC_COMM_WORLD,fileResult,"\n");
   PetscFPrintf(PETSC_COMM_WORLD,fileResult,"type: ensight \n");
   PetscFPrintf(PETSC_COMM_WORLD,fileResult,"\n");
   PetscFPrintf(PETSC_COMM_WORLD,fileResult,"GEOMETRY \n");
   PetscFPrintf(PETSC_COMM_WORLD,fileResult,"\n");
   PetscFPrintf(PETSC_COMM_WORLD,fileResult,"model: 1 Meshless.*.geo \n");
   PetscFPrintf(PETSC_COMM_WORLD,fileResult,"\n");
   PetscFPrintf(PETSC_COMM_WORLD,fileResult,"VARIABLE \n");
   PetscFPrintf(PETSC_COMM_WORLD,fileResult,"\n");
   PetscFPrintf(PETSC_COMM_WORLD,fileResult,"scalar per node: 1 Density Density.*.res \n");
   PetscFPrintf(PETSC_COMM_WORLD,fileResult,"scalar per node: 1 Pressure Pressure.*.res \n");
   PetscFPrintf(PETSC_COMM_WORLD,fileResult,"vector per node: 1 Velocity Velocity.*.res \n");
   PetscFPrintf(PETSC_COMM_WORLD,fileResult,"\n");
   PetscFPrintf(PETSC_COMM_WORLD,fileResult,"TIME \n");
   PetscFPrintf(PETSC_COMM_WORLD,fileResult,"\n");
   PetscFPrintf(PETSC_COMM_WORLD,fileResult,"time set: 1 \n");
   PetscFPrintf(PETSC_COMM_WORLD,fileResult,"number of steps: %d \n", numSteps/par->FreqResults);
   PetscFPrintf(PETSC_COMM_WORLD,fileResult,"filename start number: 0 \n");
   PetscFPrintf(PETSC_COMM_WORLD,fileResult,"filename increment: %d \n", 1);
   PetscFPrintf(PETSC_COMM_WORLD,fileResult,"time values: \n");

   PetscInt counter1=0;
   for(i=0; i< (par->finalTime/par->timeStep);i++)
   {
	   if (i % par->FreqResults ==0){
 	    PetscFPrintf(PETSC_COMM_WORLD,fileResult,"%d ", counter1);
 	    counter1=counter1 + 1;
	   }
 	    if ((i+1) % (10*par->FreqResults) == 0) PetscFPrintf(PETSC_COMM_WORLD,fileResult,"\n  ");

   }
   fclose(fileResult);
  }

  // ##################################################
  //                  Pressure File
  // ##################################################


  sprintf(filenameResults,"./Pressure.%d.res",counter);
  fileResult=fopen(filenameResults,"wt");
  if (fileResult == NULL)
      {
  	printf("Error opening file!\n");
  	exit(1);
      }
  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"Pressure \n");
    for(i=0; i<count;i++)
    {
      PetscFPrintf(PETSC_COMM_WORLD,fileResult,"%.4e  ", par->nodos[i].nodalPressure);
      if ((i+1) % 6 == 0) PetscFPrintf(PETSC_COMM_WORLD,fileResult,"\n  ");
    }
  fclose(fileResult);

  // ##################################################
  //                  Density File
  // ##################################################


  sprintf(filenameResults,"./Density.%d.res",counter);
  fileResult=fopen(filenameResults,"wt");
  if (fileResult == NULL)
      {
        printf("Error opening file!\n");
        exit(1);
      }
  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"Density \n");
    for(i=0; i<count;i++)
    {
      PetscFPrintf(PETSC_COMM_WORLD,fileResult,"%.4e ", par->nodos[i].nodalDensity);
      if ((i+1) % 6 == 0) PetscFPrintf(PETSC_COMM_WORLD,fileResult,"\n  ");
    }

  fclose(fileResult);

  // ##################################################
  //                  Velocity File
  // ##################################################

  sprintf(filenameResults,"./Velocity.%d.res",counter);
  fileResult=fopen(filenameResults,"wt");
  if (fileResult == NULL)
      {
  	printf("Error opening file!\n");
  	exit(1);
      }
  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"Velocity \n");
  for(i=0; i<count;i++)
  {
	    PetscFPrintf(PETSC_COMM_WORLD,fileResult,"%.4e   %.4e   %.4e  ", par->puntos[i].totalPhysicalVelocity[0],
	    		par->puntos[i].totalPhysicalVelocity[1],
	    		0.0);
	    if ((i+1) % 2 == 0) PetscFPrintf(PETSC_COMM_WORLD,fileResult,"\n  ");
  }
  fclose(fileResult);
  }
  PetscFunctionReturn(0);
}

typedef struct {
  PetscScalar rho,ux,uy,temp; //Density Primitive
} Field;



#undef __FUNCT__
#define __FUNCT__ "ForceY"
PetscErrorCode ForceY(IGA iga,PetscReal t,Vec U,AppCtx *user)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  DM da;
  PetscInt dof = iga->dof;
  PetscInt dim = iga->dim;

  ierr = IGACreateNodeDM(iga,dof,&da);CHKERRQ(ierr);
  Field **u;
  ierr = DMDAVecGetArray(da,U,&u);CHKERRQ(ierr);
  DMDALocalInfo info;
  ierr = DMDAGetLocalInfo(da,&info);CHKERRQ(ierr);

  PetscInt i,j;
  for(i=info.xs;i<info.xs+info.xm;i++){
    for(j=info.ys;j<info.ys+info.ym;j++){
        PetscReal x = iga->geometryX[((j-info.ys)*iga->geom_gwidth[0]+i-info.xs)*dim];
        PetscReal y = iga->geometryX[((j-info.ys)*iga->geom_gwidth[0]+i-info.xs)*dim + 1];
        u[j][i].uy       =  0.0;
    }
  }
  ierr = DMDAVecRestoreArray(da,U,&u);CHKERRQ(ierr);
  ierr = DMDestroy(&da);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "FormInitialCondition"
PetscErrorCode FormInitialCondition(IGA iga,PetscReal t,Vec U,AppCtx *user)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  DM da;
  PetscInt dof = iga->dof;
  PetscInt dim = iga->dim;

  ierr = IGACreateNodeDM(iga,dof,&da);CHKERRQ(ierr);
  Field **u;
  ierr = DMDAVecGetArray(da,U,&u);CHKERRQ(ierr);
  DMDALocalInfo info;
  ierr = DMDAGetLocalInfo(da,&info);CHKERRQ(ierr);

  PetscInt i,j;
  for(i=info.xs;i<info.xs+info.xm;i++){
    for(j=info.ys;j<info.ys+info.ym;j++){
        PetscReal x = iga->geometryX[((j-info.ys)*iga->geom_gwidth[0]+i-info.xs)*dim];
        PetscReal y = iga->geometryX[((j-info.ys)*iga->geom_gwidth[0]+i-info.xs)*dim + 1];

        PetscReal dens0  = 1659.0;
        PetscReal P0     = 100000.0;
        PetscReal A      = 495.1e9;
        PetscReal B      = 7.21e9;
        PetscReal C      = 1.62e9;
        PetscReal R1     = 4.387;
        PetscReal R2     = 0.9954;
        PetscReal omega  = 0.3469;
        PetscReal E0     = 5877.9e3;
        PetscReal Pcr    = 2.0e11;


      /*  PetscReal B = 3.31e8; //Cavitation Problem
        PetscReal P0 = 100000;
        PetscReal N = 7.15;
        PetscReal dens0 = 1000.0;

        u[j][i].ux   =  -5.55136719;
        u[j][i].uy   =  0.0;
        u[j][i].rho  =  dens0*pow((1/B)*(2.748e6-P0)+1, 1/N);
        u[j][i].temp =  290.0;

  if (x <= 3.50)
	{
    u[j][i].ux   =  -2.3820e-7;
    u[j][i].uy   =  0.0;
    u[j][i].rho  =  990.0;//1000*pow((1/B)*(22.02726-P0)+1, 1/N);
    u[j][i].temp =  290.0;
  }*/

  /*  PetscReal B = 3.31e8; //1D RDX Detonation Problem
    PetscReal P0 = 100000;
    PetscReal N = 7.15;*/
    PetscReal densWat = 1000.0;
    u[j][i].ux   =  0.0;
    u[j][i].uy   =  0.0;
    u[j][i].rho  =  densWat;
    u[j][i].temp =  290.0;

if (x <= 0.25/2+0.005103 && x >= 0.25/2-0.005103)
{
   u[j][i].ux       =  0.0;
   u[j][i].uy       =  0.0;
   u[j][i].rho      =  user->rhoRDX;
   PetscReal nu     = u[j][i].rho/dens0;
   PetscReal Ptest  = A*exp(-R1*nu) + B*exp(-R2*nu) + C/(pow(nu, 1+omega));
   PetscReal P      = Pcr;

   if(Ptest > Pcr){
     P        = Ptest;
   }else{
     P        = A*(1-omega/(R1*nu))*exp(-R1*nu) + B*(1-omega/(R2*nu))*exp(-R2*nu) + omega*dens0*E0/nu;
   }
   u[j][i].temp = P/286.6/1659.0; //Initialize temperature as ideal gas temperature
}

    }
  }
  ierr = DMDAVecRestoreArray(da,U,&u);CHKERRQ(ierr);
  ierr = DMDestroy(&da);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "OutputMonitor"
PetscErrorCode OutputMonitor(TS ts,PetscInt it_number,PetscReal c_time,Vec U,void *mctx)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  AppCtx *user = (AppCtx *)mctx;
  char           filename[256];
  sprintf(filename,"./nsvms%d.dat",it_number);
  ierr = IGAWriteVec(user->iga,U,filename);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TSSetRadius_GeneralizedAlpha"
PetscErrorCode TSSetRadius_GeneralizedAlpha(AppCtx *user,PetscReal radius)
{
  PetscFunctionBegin;
  user->Alpha_m = (3.0-radius)/(2.0*(1.0 + radius));  //(2-radius)/(1+radius);    //(3.0-radius)/(2.0*(1.0 + radius));   //
  user->Alpha_f = 1/(1+radius);
  user->Gamma   = 0.5 + user->Alpha_m - user->Alpha_f;
  user->Beta    = 0.5 * (1 + user->Alpha_m - user->Alpha_f); user->Beta *= user->Beta; //Do we need beta?
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TSPredictStage_GeneralizedAlpha"
static PetscErrorCode TSPredictStage_GeneralizedAlpha(AppCtx *user, PARAMETERS *par)
{

  Vec             V0 = user->V0, A0 = user->A0;
  Vec            Vp = user->Vp, Ap = user->Ap;
  PetscReal      dt = user->timeStep;
  PetscReal      Gamma   = user->Gamma;
  PetscReal      Beta    = user->Beta;
  PetscErrorCode ierr;
  PetscFunctionBegin;

/*
  PetscScalar *xx, *xy, *xz;
  ierr  = VecGetArray(th->X0,&xx);CHKERRQ(ierr);
  ierr  = VecGetArray(th->V0,&xy);CHKERRQ(ierr);
  ierr  = VecGetArray(th->A0,&xz);CHKERRQ(ierr);
  PetscInt k,i;
  for(i=0;i<841;i++) {
	  for(k=0;k<5;k++) {
		    PetscPrintf(PETSC_COMM_WORLD,"i: %d, X0: %e, V0: %e, A0: %e \n",i, xx[i*5+k], xy[i*5+k], xz[i*5+k]);
	  }
  }

*/

//FOR FLUID MECHANICS AND FSI
  //Vp = V0
  ierr = VecCopy(V0,Vp);CHKERRQ(ierr);     						//ok

  //Ap = (Gamma - 1)/Gamma*A0
  ierr = VecCopy(A0,Ap);CHKERRQ(ierr);
  ierr = VecScale(Ap,(Gamma - 1.0)/Gamma);CHKERRQ(ierr);		//ok

  //Xp = X0 + dt*V0 + dt²/2*((1-2*Beta)*A0 + 2*Beta*Ap)

/*
//FORN NONLINEAR SOLID MECHANICS
  // Xp = X
  ierr = VecCopy(X0,Xp);CHKERRQ(ierr);

  // Ap = - 1/(Beta*dt)*V0 - (1-2*Beta)/(2*Beta)*A0
  ierr = VecCopy(V0,Ap);CHKERRQ(ierr);  //Ap=V0
  ierr = VecAXPBY(Ap,-(1-2*Beta)/(2*Beta),-1/(dt*Beta),A0);CHKERRQ(ierr); //Ap= - 1/(Beta*dt)*Ap - (1-2*Beta)/(2*Beta)*A0


  // Vp = V0 + dt*((1-Gamma)*A0 + Gamma*Ap)
  ierr = VecWAXPY(Vp,(1.0-Gamma)/Gamma,A0,Ap);CHKERRQ(ierr); //Vp = (1.0-Gamma)/Gamma*A0 + Ap
  ierr = VecAYPX (Vp,dt*Gamma,V0);CHKERRQ(ierr); 			 //Vp = V0 + Gamma*dt*Vp
*/


  ierr = VecCopy(user->Vp,user->V1);CHKERRQ(ierr);
  ierr = VecCopy(user->Ap,user->A1);CHKERRQ(ierr);

/*
  PetscScalar *xx, *xy, *xz;
  ierr  = VecGetArray(th->X1,&xx);CHKERRQ(ierr);
  ierr  = VecGetArray(th->V1,&xy);CHKERRQ(ierr);
  ierr  = VecGetArray(th->A1,&xz);CHKERRQ(ierr);
  PetscInt k,i;
  for(i=0;i<10;i++) {
	  for(k=0;k<5;k++) {
		    PetscPrintf(PETSC_COMM_WORLD,"i: %d, Xp: %e, Vp: %e, Ap: %e \n",i, xx[i*5+k], xy[i*5+k], xz[i*5+k]);
	  }
  }
*/
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TSUpdateAlphaLevels_GeneralizedAlpha"
static PetscErrorCode TSUpdateAlphaLevels_GeneralizedAlpha(AppCtx *user, PARAMETERS *par)
{

//	PetscPrintf(PETSC_COMM_WORLD,"					TSUpdateAlphaLevels_JB\n");

  Vec            V1 = user->V1, A1 = user->A1;
  Vec            Va = user->Va, Aa = user->Aa;
  Vec            V0 = user->V0, A0 = user->A0;
  PetscReal      Alpha_m = user->Alpha_m;
  PetscReal      Alpha_f = user->Alpha_f;
  PetscErrorCode ierr;
  PetscFunctionBegin;

  /* Va = V0 + Alpha_f*(V1-V0) */
  ierr = VecWAXPY(Va,-1.0,V0,V1);CHKERRQ(ierr);
  ierr = VecAYPX (Va,Alpha_f,V0);CHKERRQ(ierr);		//ok
  /* Aa = A0 + Alpha_m*(A1-A0) */
  ierr = VecWAXPY(Aa,-1.0,A0,A1);CHKERRQ(ierr);
  ierr = VecAYPX (Aa,Alpha_m,A0);CHKERRQ(ierr);		//ok

//  PetscScalar *xx, *xy, *xz;
//  ierr  = VecGetArray(Xa,&xx);CHKERRQ(ierr);
//  ierr  = VecGetArray(Va,&xy);CHKERRQ(ierr);
//  ierr  = VecGetArray(Aa,&xz);CHKERRQ(ierr);
//  PetscInt k,i;
//  for(i=0;i<841;i++) {
//	  for(k=0;k<5;k++) {
//		    PetscPrintf(PETSC_COMM_WORLD,"i: %d, Xa_upd: %e, Va_upd: %e, Aa_upd: %e \n",i, xx[i*5+k], xy[i*5+k], xz[i*5+k]);
//	  }
//  }

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TSUpdateStage_GeneralizedAlpha"
static PetscErrorCode TSUpdateStage_GeneralizedAlpha(AppCtx *user,Vec dA)
{

//	PetscPrintf(PETSC_COMM_WORLD,"					TSUpdateStage_JB\n");

  Vec            V1 = user->V1, A1 = user->A1;
  PetscReal      dt = user->timeStep;
  PetscReal      Gamma   = user->Gamma;
  PetscReal      Beta    = user->Beta;
  PetscErrorCode ierr;

  PetscFunctionBegin;

//  PetscInt k,i;
//  PetscScalar *xx, *xy, *xz, *zz;
//  ierr  = VecGetArray(X1,&xx);CHKERRQ(ierr);
//  ierr  = VecGetArray(V1,&xy);CHKERRQ(ierr);
//  ierr  = VecGetArray(A1,&xz);CHKERRQ(ierr);
//  ierr  = VecGetArray(dA,&zz);CHKERRQ(ierr);
////   PetscInt k,i;
//  for(i=0;i<10;i++) {
//	  for(k=0;k<5;k++) {
//		    PetscPrintf(PETSC_COMM_WORLD,"i: %d, X1_b_upd: %e, V1_b_upd: %e, A1_b_upd: %e, dA_b_upd: %e \n",i, xx[i*5+k], xy[i*5+k], xz[i*5+k], zz[i*5+k]);
//	  }
//  }
//  PetscPrintf(PETSC_COMM_WORLD,"Gamma: %e, Alpha_f: %e Alpha_m: %e \n",th->Gamma, th->Alpha_f, th->Alpha_m);

//    ierr = VecDuplicate(Y,dA);
//    ierr = VecScale(dA,-1.0);CHKERRQ(ierr);
 //USING dA
    //A1 = A1 + dA
    ierr = VecAXPY (A1,-1.0,dA);CHKERRQ(ierr);  			//ok
    //V1 = V1 + Gamma*dt*dA
    ierr = VecAXPY (V1,-Gamma*dt,dA);CHKERRQ(ierr);			//ok
    //X1 = X1 + Beta*dt²*dA

    /*

//USING A
  //A1 = A

  //V1 = Vp + Gamma*dt*(A1 - A0)
  ierr = VecWAXPY(V1,-1.0,A0,A1);CHKERRQ(ierr); //V1 = (A1 - A0)
  ierr = VecAYPX (V1,Gamma*dt,Vp);CHKERRQ(ierr);
  //X1 = Xp + Beta*dt²*(A1 - A0)
  ierr = VecWAXPY(X1,-1.0,A0,A1);CHKERRQ(ierr); //X1 = (A1 - A0)
  ierr = VecAYPX (X1,Beta*dt*dt,Xp);CHKERRQ(ierr);
*/
  /*

//SOLVING X
  // A1 = Ap + 1/(dt²*Beta)*(X1-X0)
  ierr = VecWAXPY(A1,-1.0,X0,X1);CHKERRQ(ierr); //A1 = (X1 - X0)
  ierr = VecAYPX (A1,1/(dt*dt*Beta),Ap);CHKERRQ(ierr); // A1= 1/(dt*dt*Beta)*A1 + Ap

  // V1 = Vp + Gamma*dt/(dt²*Beta)*(X1 - X0)
  ierr = VecWAXPY(V1,-1.0,X0,X1);CHKERRQ(ierr); //V1 = (X1 - X0)
  ierr = VecAYPX (V1,Gamma/(dt*Beta),Vp);CHKERRQ(ierr); // V1= Gamma/(dt*Beta)*V1 + Vp

*/

//    PetscInt k,i;
//    PetscScalar *xx, *xy, *xz, *zz;
//    PetscPrintf(PETSC_COMM_WORLD,"Gamma: %e, Alpha_f: %e Alpha_m: %e \n",th->Gamma, th->Alpha_f, th->Alpha_m);

//    ierr  = VecGetArray(X1,&xx);CHKERRQ(ierr);
//    ierr  = VecGetArray(V1,&xy);CHKERRQ(ierr);
//    ierr  = VecGetArray(A1,&xz);CHKERRQ(ierr);
//    ierr  = VecGetArray(dA,&zz);CHKERRQ(ierr);
//
//    for(i=0;i<10;i++) {
//  	  for(k=0;k<5;k++) {
//  //		    PetscPrintf(PETSC_COMM_WORLD,"i: %d, X1_a_upd: %e, V1_a_upd: %e, A1_a_upd: %e, dA_a_upd: %e \n",i, xx[i*5+k], xy[i*5+k], xz[i*5+k], zz[i*5+k]);
//  	  }
//    }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeResNorm"
static PetscErrorCode ComputeResNorm(IGA iga,Vec Res, PetscInt it)
{

  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscInt		dof = iga->dof;
  PetscInt		dim = iga->dim;
  PetscReal           fnorm;
  PetscInt i,j,Nloc;
  VecGetLocalSize(Res,&Nloc);
  PetscInt Nnodes = iga->node_sizes[0]*iga->node_sizes[1]*iga->node_sizes[2];

      //Global Residual norm (including all dof)

	  ierr = VecNormBegin(Res,NORM_2,&fnorm);CHKERRQ(ierr);        // fnorm <- ||F||
	  ierr = VecNormEnd(Res,NORM_2,&fnorm);CHKERRQ(ierr);
//	  PetscPrintf(PETSC_COMM_WORLD,"it: %d, Res norm: %e  \n",it,fnorm);

	  //Residual norm by dof

	   PetscReal           normCont;
	   PetscReal           normMom;
	   PetscReal           normEner;
	   Vec                 cont;

	   VecCreate(PETSC_COMM_WORLD,&cont);
	   VecSetSizes(cont,Nloc/dof,Nnodes);
	   VecSetFromOptions(cont);

	   Vec mom;
	   VecCreate(PETSC_COMM_WORLD,&mom);
	   VecSetSizes(mom,Nloc/dof*dim,Nnodes*dim);
	   VecSetFromOptions(mom);

	   Vec ener;
	   VecCreate(PETSC_COMM_WORLD,&ener);
	   VecSetSizes(ener,Nloc/dof,Nnodes);
	   VecSetFromOptions(ener);

	   PetscScalar *arrayRes;
	   PetscScalar *aCont;
	   PetscScalar *aMom;
	   PetscScalar *aEner;

	   ierr  = VecGetArray(Res,&arrayRes);CHKERRQ(ierr);

	   ierr  = VecGetArray(cont,&aCont);CHKERRQ(ierr);
	   ierr  = VecGetArray(mom,&aMom);CHKERRQ(ierr);
	   ierr  = VecGetArray(ener,&aEner);CHKERRQ(ierr);


	      for(i=0;i<Nloc/dof;i++) {
	    	aCont[i]    = arrayRes[dof*i];
	    	aEner[i]    = arrayRes[dof*i + dof - 1];
	     	  for(j=0;j<dim;j++){
	     		aMom[dim*i + j]    = arrayRes[dof*i + 1 + j];
	     	  }
	      }

	  ierr  = VecRestoreArray(ener,&aEner);CHKERRQ(ierr);
	  ierr  = VecRestoreArray(mom,&aMom);CHKERRQ(ierr);
	  ierr  = VecRestoreArray(cont,&aCont);CHKERRQ(ierr);
	  ierr  = VecRestoreArray(Res,&arrayRes);CHKERRQ(ierr);


	   ierr = VecNormBegin(cont, NORM_2, &normCont);CHKERRQ(ierr);
	   ierr = VecNormEnd(cont, NORM_2, &normCont);CHKERRQ(ierr);

	   ierr = VecNormBegin(mom, NORM_2, &normMom);CHKERRQ(ierr);
	   ierr = VecNormEnd(mom, NORM_2, &normMom);CHKERRQ(ierr);

	   ierr = VecNormBegin(ener, NORM_2, &normEner);CHKERRQ(ierr);
	   ierr = VecNormEnd(ener, NORM_2, &normEner);CHKERRQ(ierr);

      ierr = VecDestroy(&cont);CHKERRQ(ierr);
      ierr = VecDestroy(&mom);CHKERRQ(ierr);
      ierr = VecDestroy(&ener);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char *argv[]) {

  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);

  AppCtx     user;
  PARAMETERS *par;
  ierr       = PetscMalloc1(sizeof(*par),PARAMETERS,&par);CHKERRQ(ierr);
  PetscInt   i,j,k,m,l,it;

  user.mu       = 8.9e-4;//3.137e-3;//1.983e-5;
  user.lamda    = -2.0*user.mu/3.0;
  user.kappa    = 0.6; //0.024;
  user.Cv       = 4180.0;
  user.Cp       = 4181.1;
  user.p0       = 100000;
  user.Lx       = 0.25;
  user.Ly       = 0.05000;
  user.rhoRDX   = 1659.0;

  user.F[0]=0.0;
  user.F[1]=0.0;
  user.F[2]=0.0;
  user.F[3]=0.0;

  user.currentTime  = 0.0;
  user.timeStep		  = 0.5e-10;
  user.finalTime 	  = 16e-3;
  user.stepNumber 	= 0;
  user.max_its      = 3;
  user.FreqResults  = 1000;
  user.TimeRestart  = 0;
  user.StepRestart  = 0;
  user.FreqRestarts = 2000;

  par->currentTime  = user.currentTime;
  par->timeStep		  = user.timeStep;
  par->finalTime 	  = user.finalTime;
  par->stepNumber 	= user.stepNumber;
  par->stepNumber 	= user.stepNumber;
  par->FreqResults  = user.FreqResults;

PetscInt p=2, C=PETSC_DECIDE;
ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"","PhaseFieldCrystal2D Options","IGA");CHKERRQ(ierr);
ierr = PetscOptionsInt("-StepRestart","Step of the initial solution from file",__FILE__,user.StepRestart,&user.StepRestart,NULL);CHKERRQ(ierr);
ierr = PetscOptionsReal("-TimeRestart","Time of the initial solution from file",__FILE__,user.TimeRestart,&user.TimeRestart,NULL);CHKERRQ(ierr);

ierr = PetscOptionsEnd();CHKERRQ(ierr);
if (C == PETSC_DECIDE) C = p-1;

  IGA iga;
  ierr = IGACreate(PETSC_COMM_WORLD,&iga);CHKERRQ(ierr);
  ierr = IGASetDim(iga,2);CHKERRQ(ierr);
  ierr = IGASetDof(iga,4);CHKERRQ(ierr);
  user.iga = iga;

  ierr = IGARead(iga,"./Geo/Geometry.dat");CHKERRQ(ierr);
  ierr = IGASetFromOptions(iga);CHKERRQ(ierr);
  ierr = IGASetUp(iga);CHKERRQ(ierr);
  ierr = IGAWrite(iga,"igaF.dat");CHKERRQ(ierr);
  if (par->stepNumber == 0){
    ierr = input(par, &user);CHKERRQ(ierr);
  }

PetscReal t=0;
  ierr = IGACreateVec(iga,&user.V0);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&user.A0);CHKERRQ(ierr);
  ierr = VecZeroEntries(user.V0);CHKERRQ(ierr);
  ierr = VecZeroEntries(user.A0);CHKERRQ(ierr);
  if (user.StepRestart == 0){
    ierr = FormInitialCondition(iga,t,user.V0,&user);CHKERRQ(ierr);
  } else {  //If we are restarting then we read a restart file instead of applying initial conditions
    ierr = ReadLastResults(&user, par,user.V0,user.A0,user.StepRestart);CHKERRQ(ierr);
  }
  par->currentTime  = user.currentTime;
  par->timeStep		  = user.timeStep;
  par->finalTime 	  = user.finalTime;
  user.stepNumber  = user.StepRestart;
  par->stepNumber  = user.StepRestart;
  //par->currentTime = user.TimeRestart;

  ierr = IGACreateVec(iga,&user.dA);CHKERRQ(ierr);
  ierr = VecZeroEntries(user.dA);CHKERRQ(ierr);

  // Dump Initial Solution
  char filename[256];
  sprintf(filename,"velS%d.dat",user.stepNumber);
  ierr = IGAWriteVec(user.iga,user.V0,filename);CHKERRQ(ierr);
  IGABoundary bound;

  ierr =IGASetUserIFunction(iga,Residual,&user);CHKERRQ(ierr);
  ierr =IGAGetBoundary(iga,0,0,&bound);CHKERRQ(ierr);
  ierr =IGABoundarySetUserIFunction(bound,ResidualBound,&user);CHKERRQ(ierr);
  ierr =IGAGetBoundary(iga,0,1,&bound);CHKERRQ(ierr);
  ierr =IGABoundarySetUserIFunction(bound,ResidualBound,&user);CHKERRQ(ierr);
  ierr =IGAGetBoundary(iga,1,0,&bound);CHKERRQ(ierr);
  ierr =IGABoundarySetUserIFunction(bound,ResidualBound,&user);CHKERRQ(ierr);
  ierr =IGAGetBoundary(iga,1,1,&bound);CHKERRQ(ierr);
  ierr =IGABoundarySetUserIFunction(bound,ResidualBound,&user);CHKERRQ(ierr);
  ierr =IGASetUserIJacobian(iga,Tangent,&user);CHKERRQ(ierr);

  ierr = TSSetRadius_GeneralizedAlpha(&user,0.5);CHKERRQ(ierr);

  ierr = VecDuplicate(user.A0,&user.A1);CHKERRQ(ierr);
  ierr = VecDuplicate(user.A0,&user.Ap);CHKERRQ(ierr);
  ierr = VecDuplicate(user.A0,&user.Aa);CHKERRQ(ierr);
  ierr = VecDuplicate(user.A0,&user.V1);CHKERRQ(ierr);
  ierr = VecDuplicate(user.A0,&user.Vp);CHKERRQ(ierr);
  ierr = VecDuplicate(user.A0,&user.Va);CHKERRQ(ierr);

  PetscScalar bc_p[4]={0.0};
  bc_p[0]=1.0;
  PetscScalar bc_ux[4]={0.0};
  bc_ux[1]=1.0;
  PetscScalar bc_uy[4]={0.0};
  bc_uy[2]=1.0;
  PetscScalar bc_temp[4]={0.0};
  bc_temp[3]=1.0;

  for(i=0;i<par->numNodes;i++){
    par->puntos[i].totalPhysicalVelocity[0] = 0.0;
    par->puntos[i].totalPhysicalVelocity[1] = 0.0;
    par->puntos[i].totalPhysicalAcceleration[0] = 0.0;
    par->puntos[i].totalPhysicalAcceleration[1] = 0.0;
    par->puntos[i].totalPhysicalDisplacement[0] = 0.0;
    par->puntos[i].totalPhysicalDisplacement[1] = 0.0;

  par->nodos[i].currentDeformationGradient[0] = 1.0;
  par->nodos[i].currentDeformationGradient[1] = 0.0;
  par->nodos[i].currentDeformationGradient[2] = 0.0;
  par->nodos[i].currentDeformationGradient[3] = 1.0;
  if(user.StepRestart == 0){
  par->nodos[i].nodalDensity = 0.0;
  par->nodos[i].nodalDensityInitial = 0.0;
  par->nodos[i].nodalPressure = 0.0;
  }
}

  ierr = outputTXT(par);CHKERRQ(ierr);

  const PetscScalar *arrayU;
  const PetscScalar *arrayV;
  Vec               localV;
  Vec               localU;

  //######################################
  //         	TIME STEP
  //######################################
  PetscReal Beta  = user.Beta;
  PetscReal Gamma = user.Gamma;
  PetscReal deltat= user.timeStep;
  while(user.currentTime+user.timeStep <= user.finalTime)
   {
	  	 //PetscPrintf(PETSC_COMM_WORLD,"######################################################### \n");
	  	 PetscPrintf(PETSC_COMM_WORLD,"Step Number: %d  Time step: %e, Time: %e \n",user.stepNumber, user.timeStep, user.currentTime);
	  	 //PetscPrintf(PETSC_COMM_WORLD,"######################################################### \n");
       PetscInt maxits = user.max_its;               /* maximum number of iterations */
       par->currentTime+=par->timeStep;

       //Psuedo-adaptive timestepping, enforce at least O(1) convergence
       /*if(par->stepNumber>10){
         if(Ynorm < 0.5*Ynorm0){
           user.timeStep = user.timeStep*1.01;
           par->timeStep = user.timeStep;
           deltat= user.timeStep;
        }
        if(Ynorm >= 0.5*Ynorm0 && user.timeStep > 0.1e-11){
         user.timeStep = user.timeStep*0.97;
         par->timeStep = user.timeStep;
         deltat= user.timeStep;
       }
     }*/

        //Start with initializing the particles and background for generalizedAlpha
	      ierr = TSPredictStage_GeneralizedAlpha(&user, par);CHKERRQ(ierr);

          // Nodal intermediate step
           for(i=0;i<par->numNodes;i++){
             par->puntos[i].totalPhysicalVelocityOldStep[0] = par->puntos[i].totalPhysicalVelocity[0];
             par->puntos[i].totalPhysicalVelocityOldStep[1] = par->puntos[i].totalPhysicalVelocity[1];
             par->puntos[i].totalPhysicalAccelerationOldStep[0] = par->puntos[i].totalPhysicalAcceleration[0];
             par->puntos[i].totalPhysicalAccelerationOldStep[1] = par->puntos[i].totalPhysicalAcceleration[1];
             par->puntos[i].totalPhysicalDisplacementOldStep[0] = par->puntos[i].totalPhysicalDisplacement[0];
             par->puntos[i].totalPhysicalDisplacementOldStep[1] = par->puntos[i].totalPhysicalDisplacement[1];
             par->nodos[i].DeformationGradientOld[0] = par->nodos[i].currentDeformationGradient[0];
             par->nodos[i].DeformationGradientOld[1] = par->nodos[i].currentDeformationGradient[1];
             par->nodos[i].DeformationGradientOld[2] = par->nodos[i].currentDeformationGradient[2];
             par->nodos[i].DeformationGradientOld[3] = par->nodos[i].currentDeformationGradient[3];


               for(int k = 0; k<2; k++){
                 par->puntos[i].totalPhysicalVelocity[k] =  par->puntos[i].totalPhysicalVelocityOldStep[k];
                 par->puntos[i].totalPhysicalAcceleration[k] =  ((Gamma-1)/Gamma)*par->puntos[i].totalPhysicalAccelerationOldStep[k];
                 par->puntos[i].totalPhysicalDisplacement[k] =  par->puntos[i].totalPhysicalDisplacementOldStep[k]+deltat*par->puntos[i].totalPhysicalVelocity[k]+
                 (deltat*deltat/(2.0))*((1-2*Beta)*par->puntos[i].totalPhysicalAccelerationOldStep[k]+2*Beta*par->puntos[i].totalPhysicalAcceleration[k]);
               }

             for(j=0;j<2;j++){
               par->puntos[i].tempCoord[j] = par->puntos[i].currentCoord[j];
             }
           }

	         PetscInt dof = iga->dof;
	         PetscInt dim = iga->dim;
	         Mat A0;
	         Vec Res, SP, ResFS;
           Mat MassFS;
           PetscReal nodalVolume;
           ierr = IGACreateMat(iga,&MassFS);CHKERRQ(ierr);
	         ierr = IGACreateMat(iga,&A0);CHKERRQ(ierr);
	         ierr = IGACreateVec(iga,&Res);CHKERRQ(ierr);
	         ierr = IGACreateVec(iga,&SP);CHKERRQ(ierr);
           ierr = IGACreateVec(iga,&ResFS);CHKERRQ(ierr);

	         for (it=0; it<maxits; it++) {
	       	  //PetscPrintf(PETSC_COMM_WORLD,":::::::::::::::::::::::::::::::::\n");
	       	  //PetscPrintf(PETSC_COMM_WORLD,"  Iteration: %d  \n", it);
	       	  //PetscPrintf(PETSC_COMM_WORLD,":::::::::::::::::::::::::::::::::\n");

	       	  ierr = TSUpdateAlphaLevels_GeneralizedAlpha(&user, &par);CHKERRQ(ierr);

	             PetscReal t = user.currentTime;
	             PetscReal dt = user.timeStep;
	             PetscReal Alpha_m = user.Alpha_m;
	             PetscReal Alpha_f = user.Alpha_f;
               par->Alpha_f = Alpha_f;
	             PetscReal stage_time = t + Alpha_f*dt;

           for(i=0;i<par->numNodes;i++){
               par->puntos[i].totalPhysicalVelocityOldIteration[0] = par->puntos[i].totalPhysicalVelocity[0];
               par->puntos[i].totalPhysicalVelocityOldIteration[1] = par->puntos[i].totalPhysicalVelocity[1];
               par->puntos[i].totalPhysicalAccelerationOldIteration[0] = par->puntos[i].totalPhysicalAcceleration[0];
               par->puntos[i].totalPhysicalAccelerationOldIteration[1] = par->puntos[i].totalPhysicalAcceleration[1];
               par->puntos[i].totalPhysicalDisplacementOldIteration[0] = par->puntos[i].totalPhysicalDisplacement[0];
               par->puntos[i].totalPhysicalDisplacementOldIteration[1] = par->puntos[i].totalPhysicalDisplacement[1];
               par->puntos[i].totalPhysicalVelocity[0] = par->puntos[i].totalPhysicalVelocityOldStep[0] + Alpha_f*(par->puntos[i].totalPhysicalVelocity[0] - par->puntos[i].totalPhysicalVelocityOldStep[0]);
               par->puntos[i].totalPhysicalDisplacement[0] = par->puntos[i].totalPhysicalDisplacementOldStep[0] + Alpha_f*(par->puntos[i].totalPhysicalDisplacement[0] - par->puntos[i].totalPhysicalDisplacementOldStep[0]);
               par->puntos[i].totalPhysicalVelocity[1] = par->puntos[i].totalPhysicalVelocityOldStep[1] + Alpha_f*(par->puntos[i].totalPhysicalVelocity[1] - par->puntos[i].totalPhysicalVelocityOldStep[1]);
               par->puntos[i].totalPhysicalDisplacement[1] = par->puntos[i].totalPhysicalDisplacementOldStep[1] + Alpha_f*(par->puntos[i].totalPhysicalDisplacement[1] - par->puntos[i].totalPhysicalDisplacementOldStep[1]);
             for(j=0;j<2;j++){
                par->puntos[i].currentCoord[j] = par->puntos[i].tempCoord[j] + par->puntos[i].totalPhysicalDisplacement[j] - par->puntos[i].totalPhysicalDisplacementOldStep[j];
              }
              }


              Vec Mass;
              ierr = IGACreateVec(iga,&Mass);CHKERRQ(ierr);
              ierr = VecZeroEntries(Mass);CHKERRQ(ierr);

              PetscReal factor = 1.0;
              ierr = IGAGetLocalVecArray(iga,user.Va,&localU,&arrayU);CHKERRQ(ierr);
              ierr = IGAGetLocalVecArray(iga,user.Aa,&localV,&arrayV);CHKERRQ(ierr);

              ierr = ComputeCurrentExplosiveVolume_dens(&user, par, &arrayV[0],&arrayU[0], it);CHKERRQ(ierr);
              //PetscPrintf(PETSC_COMM_WORLD,"Explosive Volume Change =  %e \n", user.totalCurrentExplosiveVolume - user.totalInitialExplosiveVolume);

             ierr = VecZeroEntries(Res);CHKERRQ(ierr);
	           ierr = IGAComputeIFunction(iga,dt,Alpha_m,user.Aa,stage_time,user.Va,Res);CHKERRQ(ierr);

             ierr = VecZeroEntries(ResFS);CHKERRQ(ierr);
             ierr = IGAComputeResidualFS(par,iga,dt,Alpha_m,&arrayV[0],stage_time,&arrayU[0],ResFS,&user);CHKERRQ(ierr);
             ierr = VecAssemblyBegin(ResFS);CHKERRQ(ierr);
             ierr = VecAssemblyEnd(ResFS);CHKERRQ(ierr);

             ierr = VecAXPY(Res,factor,ResFS);CHKERRQ(ierr);
             ierr = MatZeroEntries(MassFS);CHKERRQ(ierr);
             ierr = IGAComputeIJacobianFS(Alpha_m,MassFS,&user,&arrayU[0], par);CHKERRQ(ierr);
	           ierr = IGAComputeIJacobianComp(iga,dt,Alpha_m,user.Aa,stage_time,user.Va,A0,&user);CHKERRQ(ierr); //Integral over the whole domain of the
             //Immersing fluid, from which the MassFS can be subtracted to get only the fluid component plus the solid component
             ierr = MatAXPY(A0,-factor,MassFS,SAME_NONZERO_PATTERN);CHKERRQ(ierr); //subtract integral of N_A*[A^f-A^RDX]*N_B
             // on foreground

	            PetscInt numFluidNodes  = iga->geom_lwidth[0]*iga->geom_lwidth[1]*iga->geom_lwidth[2];
    	        PetscInt nodesX  = iga->geom_lwidth[0], nodesY  = iga->geom_lwidth[1], nodesZ  = iga->geom_lwidth[2];
    	        PetscInt gnodesX = iga->geom_gwidth[0], gnodesY = iga->geom_gwidth[1];

      	         for(m=0;m<nodesZ;m++) {
      	             for(l=0;l<nodesY;l++) {
      	                 for(k=0;k<nodesX;k++) {
      	            /*if(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim + 0] <=0.00001) {
      	            	PetscInt index = (m*gnodesX*gnodesY + l*gnodesX+ k)*dof;
      	            	PetscInt index_array[4]={0.0};
      	            	index_array[0]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof;
      	            	index_array[1]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+1;
      	            	index_array[2]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+2;
      	            	index_array[3]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+3;

      	            	MatSetValuesLocal(A0,1,&index,4,index_array,bc_p,INSERT_VALUES);CHKERRQ(ierr);
      	            	MatSetValuesLocal(A0,4,index_array,1,&index,bc_p,INSERT_VALUES);CHKERRQ(ierr);
      	            	VecSetValueLocal(Res,index,0.0,INSERT_VALUES);CHKERRQ(ierr);

      	            }*/
                    if(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim + 0] <=0.00001) {
                      PetscInt index = (m*gnodesX*gnodesY + l*gnodesX+ k)*dof+1;
                      PetscInt index_array[4]={0.0};
                      index_array[0]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof;
                      index_array[1]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+1;
                      index_array[2]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+2;
                      index_array[3]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+3;
                      MatSetValuesLocal(A0,1,&index,4,index_array,bc_ux,INSERT_VALUES);CHKERRQ(ierr);
                      MatSetValuesLocal(A0,4,index_array,1,&index,bc_ux,INSERT_VALUES);CHKERRQ(ierr);
                      VecSetValueLocal(Res,index,0.0,INSERT_VALUES);CHKERRQ(ierr);

                    }

     	            if(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim] >=(user.Lx -0.00001)) {
      	             	            	PetscInt index = (m*gnodesX*gnodesY + l*gnodesX+ k)*dof+1;
      	             	            	PetscInt index_array[4]={0.0};
      	             	            	index_array[0]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof;
      	             	            	index_array[1]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+1;
      	             	            	index_array[2]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+2;
      	             	            	index_array[3]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+3;
      	             	            	MatSetValuesLocal(A0,1,&index,4,index_array,bc_ux,INSERT_VALUES);CHKERRQ(ierr);
      	             	            	MatSetValuesLocal(A0,4,index_array,1,&index,bc_ux,INSERT_VALUES);CHKERRQ(ierr);
      	             	            	VecSetValueLocal(Res,index,0.0,INSERT_VALUES);CHKERRQ(ierr);

      	             	            }
                  /*if(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim] >=(user.Lx -0.00001)) {
                      	            PetscInt index = (m*gnodesX*gnodesY + l*gnodesX+ k)*dof;
                      	            PetscInt index_array[4]={0.0};
                      	            index_array[0]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof;
                      	            index_array[1]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+1;
                      	            index_array[2]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+2;
                      	            index_array[3]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+3;

                      	            MatSetValuesLocal(A0,1,&index,4,index_array,bc_p,INSERT_VALUES);CHKERRQ(ierr);
                      	            MatSetValuesLocal(A0,4,index_array,1,&index,bc_p,INSERT_VALUES);CHKERRQ(ierr);
                      	            VecSetValueLocal(Res,index,0.0,INSERT_VALUES);CHKERRQ(ierr);
                  }*/

                  if(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim] >=(user.Lx -0.00001)) {
                      	            PetscInt index = (m*gnodesX*gnodesY + l*gnodesX+ k)*dof+3;
                      	            PetscInt index_array[4]={0.0};
                      	            index_array[0]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof;
                      	            index_array[1]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+1;
                      	            index_array[2]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+2;
                      	            index_array[3]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+3;
                      	            MatSetValuesLocal(A0,1,&index,4,index_array,bc_temp,INSERT_VALUES);CHKERRQ(ierr);
                      	            MatSetValuesLocal(A0,4,index_array,1,&index,bc_temp,INSERT_VALUES);CHKERRQ(ierr);
                      	            VecSetValueLocal(Res,index,0.0,INSERT_VALUES);CHKERRQ(ierr);
                  }
                  if(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim] <= (0.00001)) {
                      	            PetscInt index = (m*gnodesX*gnodesY + l*gnodesX+ k)*dof+3;
                      	            PetscInt index_array[4]={0.0};
                      	            index_array[0]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof;
                      	            index_array[1]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+1;
                      	            index_array[2]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+2;
                      	            index_array[3]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+3;
                      	            MatSetValuesLocal(A0,1,&index,4,index_array,bc_temp,INSERT_VALUES);CHKERRQ(ierr);
                      	            MatSetValuesLocal(A0,4,index_array,1,&index,bc_temp,INSERT_VALUES);CHKERRQ(ierr);
                      	            VecSetValueLocal(Res,index,0.0,INSERT_VALUES);CHKERRQ(ierr);
                  }

                  if(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim +1] <=0.00001) {
                    PetscInt index = (m*gnodesX*gnodesY + l*gnodesX+ k)*dof+2;
                    PetscInt index_array[4]={0.0};
                    index_array[0]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof;
                    index_array[1]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+1;
                    index_array[2]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+2;
                    index_array[3]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+3;
                    MatSetValuesLocal(A0,1,&index,4,index_array,bc_uy,INSERT_VALUES);CHKERRQ(ierr);
                    MatSetValuesLocal(A0,4,index_array,1,&index,bc_uy,INSERT_VALUES);CHKERRQ(ierr);
                    VecSetValueLocal(Res,index,0.0,INSERT_VALUES);CHKERRQ(ierr);
                  }

                  if(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim +1] >=(user.Ly -0.00001)) {
                    PetscInt index = (m*gnodesX*nodesY + l*gnodesX+ k)*dof+2;
                    PetscInt index_array[4]={0.0};
                    index_array[0]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof;
                    index_array[1]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+1;
                    index_array[2]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+2;
                    index_array[3]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+3;
                    MatSetValuesLocal(A0,1,&index,4,index_array,bc_uy,INSERT_VALUES);CHKERRQ(ierr);
                    MatSetValuesLocal(A0,4,index_array,1,&index,bc_uy,INSERT_VALUES);CHKERRQ(ierr);
                    VecSetValueLocal(Res,index,0.0,INSERT_VALUES);CHKERRQ(ierr);}
      	                 }
      	             }
      	         }

    	         ierr = MatAssemblyBegin(A0,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    	         ierr = MatAssemblyEnd(A0,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    	         ierr = VecAssemblyBegin(Res);CHKERRQ(ierr);
    	         ierr = VecAssemblyEnd(Res);CHKERRQ(ierr);

    	            KSP ksp;
    	               ierr = IGACreateKSP(iga,&ksp);CHKERRQ(ierr);
    	               ierr = KSPSetOperators(ksp,A0,A0,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
    	               ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
    	               ierr = KSPSolve(ksp,Res,user.dA);CHKERRQ(ierr);

                  PetscReal Ynorm;
		                 ierr = VecNormBegin(user.dA,NORM_2,&Ynorm);CHKERRQ(ierr);
		                 ierr = VecNormEnd(user.dA,NORM_2,&Ynorm);CHKERRQ(ierr);
		                 PetscPrintf(PETSC_COMM_WORLD,"it: %d, Y norm: %e  \n",it,Ynorm);

	           ierr = TSUpdateStage_GeneralizedAlpha(&user,user.dA);CHKERRQ(ierr);
	           ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
             ierr = VecDestroy(&Mass);CHKERRQ(ierr);

             const PetscScalar *arrayU1;
             Vec               localU1;

             ierr = IGAGetLocalVecArray(iga,user.V1,&localU1,&arrayU1);CHKERRQ(ierr);
             for (i=0;i<par->numNodes;i++){
               ierr = interpolateVelocityOnSolidNodes(par,&user,i, user.V1,&arrayU1[0]);CHKERRQ(ierr);  //Projecting the velocities from the background grid to the particles.
             }
             ierr = IGARestoreLocalVecArray(iga,user.V1,&localU1,&arrayU1);CHKERRQ(ierr);
             for (i=0;i<par->numNodes;i++){
               par->puntos[i].AccelerationIncrement[0] = par->puntos[i].totalPhysicalVelocity[0] - par->puntos[i].totalPhysicalVelocityOldIteration[0];
               par->puntos[i].AccelerationIncrement[0] /= user.Gamma*par->timeStep;
               par->puntos[i].AccelerationIncrement[1] = par->puntos[i].totalPhysicalVelocity[1] - par->puntos[i].totalPhysicalVelocityOldIteration[1];
               par->puntos[i].AccelerationIncrement[1] /= user.Gamma*par->timeStep;
             }
             //Obtain solid points' acceleration at current step and iteration
             for (i=0;i<par->numNodes;i++){
               par->puntos[i].totalPhysicalAcceleration[0] = par->puntos[i].totalPhysicalAccelerationOldIteration[0] + par->puntos[i].AccelerationIncrement[0];
               par->puntos[i].totalPhysicalAcceleration[1] = par->puntos[i].totalPhysicalAccelerationOldIteration[1] + par->puntos[i].AccelerationIncrement[1];
             }
              //Update solid displacement consistent with Newmark-Beta
              for(i=0;i<par->numNodes;i++){
                 par->puntos[i].totalPhysicalDisplacement[0] = par->puntos[i].totalPhysicalDisplacementOldIteration[0];
                 par->puntos[i].totalPhysicalDisplacement[0] += user.Beta*par->timeStep*par->timeStep*par->puntos[i].AccelerationIncrement[0];

                 par->puntos[i].totalPhysicalDisplacement[1] = par->puntos[i].totalPhysicalDisplacementOldIteration[1];
                 par->puntos[i].totalPhysicalDisplacement[1] += user.Beta*par->timeStep*par->timeStep*par->puntos[i].AccelerationIncrement[1];

               }

               ierr = VecDestroy(&localU1);CHKERRQ(ierr);
               ierr = PetscFree(arrayU1);CHKERRQ(ierr);

           }  //End iteration loop

           for(i=0;i<par->numNodes;i++) {
              for(j=0;j<2;j++){
                 par->puntos[i].currentCoord[j]=par->puntos[i].tempCoord[j] + par->puntos[i].totalPhysicalDisplacement[j] - par->puntos[i].totalPhysicalDisplacementOldStep[j];
              }
           }
             ierr = ForceY(iga,t,user.V0,&user);CHKERRQ(ierr);
	           ierr = VecCopy(user.V1,user.V0);CHKERRQ(ierr);
	           ierr = VecCopy(user.A1,user.A0);CHKERRQ(ierr);

	      // Dump solution vector
	      if (user.stepNumber % user.FreqResults == 0) {
	      char filename[256];
	      sprintf(filename,"velS%d.dat",user.stepNumber);
	      ierr = IGAWriteVec(user.iga,user.V1,filename);CHKERRQ(ierr);
	      }
	      // Dump Restart vector
	      if (user.stepNumber % user.FreqRestarts == 0) {
	      char filename[256];
	      sprintf(filename,"velS%d.dat",user.stepNumber);
	      ierr = IGAWriteVec(user.iga,user.V1,filename);CHKERRQ(ierr);
	      sprintf(filename,"acelS%d.dat",user.stepNumber);
	      ierr = IGAWriteVec(user.iga,user.A1,filename);CHKERRQ(ierr);
	      }
        if (par->stepNumber % user.FreqRestarts == 0){
          ierr = OutputRestarts(par,user.V1,user.A1);CHKERRQ(ierr);
        }
        if (par->stepNumber % user.FreqRestarts == user.FreqRestarts-1){
          ierr = OutputOldGeometry(par);CHKERRQ(ierr);
        }

        user.currentTime+=user.timeStep;
  	    user.stepNumber++;
        par->stepNumber++;
        ierr = outputTXT(par);CHKERRQ(ierr);

    ierr = MatDestroy(&A0);CHKERRQ(ierr);
    ierr = VecDestroy(&Res);CHKERRQ(ierr);
    ierr = VecDestroy(&ResFS);CHKERRQ(ierr);
    ierr = MatDestroy(&MassFS);CHKERRQ(ierr);
    ierr = VecDestroy(&SP);CHKERRQ(ierr);
   }   //End time loop

    ierr = VecDestroy(&user.V0);CHKERRQ(ierr);
    ierr = VecDestroy(&user.Va);CHKERRQ(ierr);
    ierr = VecDestroy(&user.V1);CHKERRQ(ierr);

    ierr = VecDestroy(&user.A0);CHKERRQ(ierr);
    ierr = VecDestroy(&user.Aa);CHKERRQ(ierr);
    ierr = VecDestroy(&user.A1);CHKERRQ(ierr);

    ierr = VecDestroy(&user.Vp);CHKERRQ(ierr);
    ierr = VecDestroy(&user.Ap);CHKERRQ(ierr);

    ierr = VecDestroy(&user.dA);CHKERRQ(ierr);

    ierr = PetscFree(par);CHKERRQ(ierr);
    ierr = IGADestroy(&iga);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD," DONE\n");

    ierr = PetscFinalize();CHKERRQ(ierr);

  return 0;
}
