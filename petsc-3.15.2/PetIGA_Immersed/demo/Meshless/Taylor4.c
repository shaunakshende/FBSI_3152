//Taylor Bar with immersed method. Explicit generalized alpha in accelerations

#include "petiga.h"

typedef struct {
  PetscReal nu;
  PetscReal rho;
  PetscReal lambda;
  PetscReal mu;
  PetscReal E;
  PetscReal dt;
  PetscReal final_time;
  PetscInt freqResults;
  PetscInt step;
  PetscScalar *Stress; //History variable stress
  PetscScalar *EpsilonPlastEq; //History variable equivalent plastic strain
  PetscScalar *Stress0;
  PetscReal stressY; //initial yield stress
  PetscReal h; //hardening
  PetscInt it; //current iteration
  PetscInt maxIt; //Maximum number of iterations
  PetscReal Alpha_m;
  PetscReal Alpha_f;
  PetscReal Gamma;
  PetscReal Beta;
  PetscReal Lx;
  PetscReal Ly;
  PetscReal Lz;
  Vec dV;
  Vec U0;
  Vec Up;
  Vec U1;
  Vec Ua;
  Vec V0;
  Vec Vp;
  Vec V1;
  Vec Va;
  IGA iga;

  PetscInt numParticles;
  PetscInt nen;


} AppCtx;

typedef struct
{
    PetscInt  ID;
    PetscReal initialCoord[3];
    PetscReal currentCoord[3];
    PetscReal Displacement[3];
    PetscReal DisplacementOldIteration[3];
    PetscReal DisplacementOldStep[3];
    PetscReal Velocity[3];
    PetscReal VelocityOldIteration[3];
    PetscReal VelocityOldStep[3];
    PetscReal Acceleration[3];
    PetscReal AccelerationOldIteration[3];
    PetscReal AccelerationOldStep[3];
    PetscReal AccelerationIncrement[3];
    PetscReal Stress[6];
    PetscReal Stress0[6];
    PetscReal EpsilonPlastEq;
    PetscReal DeltaEpsilonPlastEq;
    PetscReal nodalVolume;
    PetscReal velocityGradient[9];

    PetscScalar *N0;
    PetscScalar *Nx;
    PetscScalar *Ny;
    PetscScalar *Nz;
    PetscInt *map;

} PARTICLES;



PetscErrorCode ReadNumParticles(void *ctx)
{
	PetscErrorCode ierr;
	AppCtx *user = (AppCtx *)ctx;

	char filenameCoord[256];
	sprintf(filenameCoord,"./input_coor1.dat");
    FILE *fL1 = fopen(filenameCoord, "rt");
    ierr = fscanf(fL1, "%d", &user->numParticles);
    fclose(fL1);

	return 0;
}

PetscErrorCode Input(void *ctx,PARTICLES *particles)
{
	PetscErrorCode ierr;
	AppCtx *user = (AppCtx *)ctx;
	PetscInt i,j;

	char filenameCoord[256];
	sprintf(filenameCoord,"./input_coor1.dat");
    FILE *fL1 = fopen(filenameCoord, "rt");
    ierr = fscanf(fL1, "%d", &user->numParticles);

    for (i=0;i<user->numParticles;i++){
    	ierr = fscanf(fL1, "%d",&particles[i].ID);
    	ierr = fscanf(fL1, "%le ",&particles[i].initialCoord[0]);
    	ierr = fscanf(fL1, "%le ",&particles[i].initialCoord[1]);
    	ierr = fscanf(fL1, "%le ",&particles[i].initialCoord[2]);
    }
    fclose(fL1);

    for (i=0;i<user->numParticles;i++){
    	particles[i].initialCoord[0] += 0.00509;
    	particles[i].initialCoord[1] += 0.00509;
    	particles[i].initialCoord[2] += 0.00001;

    	particles[i].currentCoord[0] = particles[i].initialCoord[0];
    	particles[i].currentCoord[1] = particles[i].initialCoord[1];
    	particles[i].currentCoord[2] = particles[i].initialCoord[2];
    }

	char filenameXVOL[256];
	sprintf(filenameXVOL,"./XVOL1.dat");
	FILE *fL2 = fopen(filenameXVOL, "rt");

	for (i=0;i<user->numParticles;i++){
		ierr = fscanf(fL2, "%le", &particles[i].nodalVolume);
	}
	fclose(fL2);

	for (i=0;i<user->numParticles;i++){
		particles[i].EpsilonPlastEq = 0.0;
		particles[i].DeltaEpsilonPlastEq = 0.0;
		for (j=0;j<3;j++){
			particles[i].Displacement[j] = 0.0;
			particles[i].Velocity[j] = 0.0;
			particles[i].Acceleration[j] = 0.0;
		}
		for (j=0;j<6;j++){
			particles[i].Stress[j] = 0.0;
			particles[i].Stress0[j] = 0.0;
		}
	}

	Vec localU;
	const PetscScalar *arrayU;
	ierr = IGAGetLocalVecArray(user->iga,user->U0,&localU,&arrayU);CHKERRQ(ierr);
	PetscInt n;
	ierr = VecGetSize(localU,&n);
	ierr = IGARestoreLocalVecArray(user->iga,user->U0,&localU,&arrayU);CHKERRQ(ierr);

	user->Lx = user->iga->geometryX[n-3];
	user->Ly = user->iga->geometryX[n-2];
	user->Lz = user->iga->geometryX[n-1];

	user->nen = 0;
	user->nen = user->iga->axis[0]->p + 1;
	user->nen *= user->iga->axis[1]->p + 1;
	user->nen *= user->iga->axis[2]->p + 1;
	for (i=0;i<user->numParticles;i++){
		ierr = PetscMalloc1(user->nen,&particles[i].N0);CHKERRQ(ierr);
		ierr = PetscMalloc1(user->nen,&particles[i].Nx);CHKERRQ(ierr);
		ierr = PetscMalloc1(user->nen,&particles[i].Ny);CHKERRQ(ierr);
		ierr = PetscMalloc1(user->nen,&particles[i].Nz);CHKERRQ(ierr);
		ierr = PetscMalloc1(user->nen,&particles[i].map);CHKERRQ(ierr);
	}

	return 0;
}


PetscErrorCode FormInitialCondition(AppCtx *user,IGA iga,PetscReal t,Vec U)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = VecZeroEntries(U);CHKERRQ(ierr);
  DM da;
  ierr = IGACreateNodeDM(iga,3,&da);CHKERRQ(ierr);
  DMDALocalInfo info;
  ierr = DMDAGetLocalInfo(da,&info);CHKERRQ(ierr);
  PetscScalar ****u;
  ierr = DMDAVecGetArrayDOF(da,U,&u);CHKERRQ(ierr);

  Vec localU;
  const PetscScalar *arrayU;
  ierr = IGAGetLocalVecArray(user->iga,U,&localU,&arrayU);CHKERRQ(ierr);
  PetscInt N;
  ierr = VecGetSize(localU,&N);
  ierr = IGARestoreLocalVecArray(user->iga,U,&localU,&arrayU);CHKERRQ(ierr);

  PetscReal Lz = iga->geometryX[N-1];

  PetscInt  i,j,k;
  for(k=info.zs;k<info.zs+info.zm;k++){
    for(j=info.ys;j<info.ys+info.ym;j++){
      for(i=info.xs;i<info.xs+info.xm;i++){
    	  PetscReal z = (PetscReal)k / ( (PetscReal)(info.mz-1) )*Lz;
//    	  PetscReal z = user->iga->geometryX[((k-info.zs)*iga->geom_gwidth[0]*iga->geom_gwidth[1]
//		  +(j-info.ys)*iga->geom_gwidth[0]+i-info.xs)*dim + 2];

    	  if (z<=0.000001){
    		  u[k][j][i][2] = 0.0;
    	  }else{
    		  u[k][j][i][2] = -250.0;
    	  }

      }
    }
  }

  ierr = DMDAVecRestoreArrayDOF(da,U,&u);CHKERRQ(ierr);
  ierr = DMDestroy(&da);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


PetscErrorCode SetParticlesOldStep(void *ctx,PARTICLES *particles)
{
	AppCtx *user = (AppCtx *)ctx;
	PetscInt i;

	for (i=0;i<user->numParticles;i++){
		particles[i].VelocityOldStep[0] = particles[i].Velocity[0];
		particles[i].VelocityOldStep[1] = particles[i].Velocity[1];
		particles[i].VelocityOldStep[2] = particles[i].Velocity[2];

		particles[i].AccelerationOldStep[0] = particles[i].Acceleration[0];
		particles[i].AccelerationOldStep[1] = particles[i].Acceleration[1];
		particles[i].AccelerationOldStep[2] = particles[i].Acceleration[2];

		particles[i].DisplacementOldStep[0] = particles[i].Displacement[0];
		particles[i].DisplacementOldStep[1] = particles[i].Displacement[1];
		particles[i].DisplacementOldStep[2] = particles[i].Displacement[2];
	}

	return 0;
}

PetscErrorCode PredictorStage(void *ctx,PARTICLES *particles)
{
  PetscErrorCode ierr;
  AppCtx *user = (AppCtx *)ctx;
  PetscReal dt = user->dt;
  PetscInt i;

  ierr = VecCopy(user->U0,user->Up);CHKERRQ(ierr); //Same velocity predictor

  ierr = VecCopy(user->V0,user->Vp);CHKERRQ(ierr);
  ierr = VecScale(user->Vp,(user->Gamma - 1.0)/user->Gamma);CHKERRQ(ierr);

  ierr = VecCopy(user->Up,user->U1);CHKERRQ(ierr);
  ierr = VecCopy(user->Vp,user->V1);CHKERRQ(ierr);

  PetscReal factor1 = (user->Gamma - 1) / user->Gamma;
  PetscReal factor2 = dt * dt / 2 * (1 - 2 * user->Beta);

  for (i=0;i<user->numParticles;i++){
	  particles[i].Acceleration[0] = factor1 * particles[i].Acceleration[0];
	  particles[i].Acceleration[1] = factor1 * particles[i].Acceleration[1];
	  particles[i].Acceleration[2] = factor1 * particles[i].Acceleration[2];

	  particles[i].Displacement[0] = particles[i].DisplacementOldStep[0];
	  particles[i].Displacement[0] += dt * particles[i].VelocityOldStep[0];
	  particles[i].Displacement[0] += factor2 * particles[i].AccelerationOldStep[0];
	  particles[i].Displacement[0] += dt * dt * user->Beta * particles[i].Acceleration[0];

	  particles[i].Displacement[1] = particles[i].DisplacementOldStep[1];
	  particles[i].Displacement[1] += dt * particles[i].VelocityOldStep[1];
	  particles[i].Displacement[1] += factor2 * particles[i].AccelerationOldStep[1];
	  particles[i].Displacement[1] += dt * dt * user->Beta * particles[i].Acceleration[1];

	  particles[i].Displacement[2] = particles[i].DisplacementOldStep[2];
	  particles[i].Displacement[2] += dt * particles[i].VelocityOldStep[2];
	  particles[i].Displacement[2] += factor2 * particles[i].AccelerationOldStep[2];
	  particles[i].Displacement[2] += dt * dt * user->Beta * particles[i].Acceleration[2];

  }
  return 0;
}


PetscErrorCode SetParticlesOldIteration(void *ctx,PARTICLES *particles)
{
	AppCtx *user = (AppCtx *)ctx;
	PetscInt i;

	for (i=0;i<user->numParticles;i++){
		particles[i].VelocityOldIteration[0] = particles[i].Velocity[0];
		particles[i].VelocityOldIteration[1] = particles[i].Velocity[1];
		particles[i].VelocityOldIteration[2] = particles[i].Velocity[2];

		particles[i].AccelerationOldIteration[0] = particles[i].Acceleration[0];
		particles[i].AccelerationOldIteration[1] = particles[i].Acceleration[1];
		particles[i].AccelerationOldIteration[2] = particles[i].Acceleration[2];

		particles[i].DisplacementOldIteration[0] = particles[i].Displacement[0];
		particles[i].DisplacementOldIteration[1] = particles[i].Displacement[1];
		particles[i].DisplacementOldIteration[2] = particles[i].Displacement[2];
	}

	return 0;
}


PetscErrorCode AlphaStage(void *ctx,PARTICLES *particles)
{
  PetscErrorCode ierr;
  AppCtx *user = (AppCtx *)ctx;
  PetscInt i;

  ierr = VecWAXPY(user->Va,-1.0,user->V0,user->V1);CHKERRQ(ierr);
  ierr = VecAYPX(user->Va,user->Alpha_m,user->V0);CHKERRQ(ierr);

  ierr = VecWAXPY(user->Ua,-1.0,user->U0,user->U1);CHKERRQ(ierr);
  ierr = VecAYPX(user->Ua,user->Alpha_f,user->U0);CHKERRQ(ierr);


  for (i=0;i<user->numParticles;i++){

	  particles[i].Acceleration[0] = particles[i].AccelerationOldStep[0] + user->Alpha_m * particles[i].Acceleration[0];
	  particles[i].Acceleration[0] -= user->Alpha_m * particles[i].AccelerationOldStep[0];

	  particles[i].Acceleration[1] = particles[i].AccelerationOldStep[1] + user->Alpha_m * particles[i].Acceleration[1];
	  particles[i].Acceleration[1] -= user->Alpha_m * particles[i].AccelerationOldStep[1];

	  particles[i].Acceleration[2] = particles[i].AccelerationOldStep[2] + user->Alpha_m * particles[i].Acceleration[2];
	  particles[i].Acceleration[2] -= user->Alpha_m * particles[i].AccelerationOldStep[2];

	  particles[i].Velocity[0] = particles[i].VelocityOldStep[0] + user->Alpha_f * particles[i].Velocity[0];
	  particles[i].Velocity[0] -= user->Alpha_f * particles[i].VelocityOldStep[0];

	  particles[i].Velocity[1] = particles[i].VelocityOldStep[1] + user->Alpha_f * particles[i].Velocity[1];
	  particles[i].Velocity[1] -= user->Alpha_f * particles[i].VelocityOldStep[1];

	  particles[i].Velocity[2] = particles[i].VelocityOldStep[2] + user->Alpha_f * particles[i].Velocity[2];
	  particles[i].Velocity[2] -= user->Alpha_f * particles[i].VelocityOldStep[2];

	  particles[i].Displacement[0] = particles[i].DisplacementOldStep[0] + user->Alpha_f * particles[i].Displacement[0];
	  particles[i].Displacement[0] -= user->Alpha_f * particles[i].DisplacementOldStep[0];

	  particles[i].Displacement[1] = particles[i].DisplacementOldStep[1] + user->Alpha_f * particles[i].Displacement[1];
	  particles[i].Displacement[1] -= user->Alpha_f * particles[i].DisplacementOldStep[1];

	  particles[i].Displacement[2] = particles[i].DisplacementOldStep[2] + user->Alpha_f * particles[i].Displacement[2];
	  particles[i].Displacement[2] -= user->Alpha_f * particles[i].DisplacementOldStep[2];

  }

  for (i=0;i<user->numParticles;i++){
	  particles[i].currentCoord[0] = particles[i].initialCoord[0] + particles[i].Displacement[0];
	  particles[i].currentCoord[1] = particles[i].initialCoord[1] + particles[i].Displacement[1];
	  particles[i].currentCoord[2] = particles[i].initialCoord[2] + particles[i].Displacement[2];
  }

  return 0;
}




PetscErrorCode ComputeVelocityGradientAndShapeFunctions(void *ctx,PARTICLES *particles)
{
  PetscErrorCode ierr;
  AppCtx *user = (AppCtx *)ctx;
  PetscInt i,j;
  PetscReal pt[3] = {0.0};
  PetscScalar grad_u[3][3];

  for (i=0;i<user->numParticles;i++){
	  pt[0] = particles[i].currentCoord[0] / user->Lx;
	  pt[1] = particles[i].currentCoord[1] / user->Ly;
	  pt[2] = particles[i].currentCoord[2] / user->Lz;

	  IGAProbe prb;
	  ierr = IGAProbeCreate(user->iga,user->Ua,&prb);CHKERRQ(ierr);
	  ierr = IGAProbeSetPoint(prb,pt);CHKERRQ(ierr);
	  ierr = IGAProbeFormGrad(prb,&grad_u[0][0]);CHKERRQ(ierr);


	  particles[i].velocityGradient[0] = grad_u[0][0];
	  particles[i].velocityGradient[1] = grad_u[0][1];
	  particles[i].velocityGradient[2] = grad_u[0][2];

	  particles[i].velocityGradient[3] = grad_u[1][0];
	  particles[i].velocityGradient[4] = grad_u[1][1];
	  particles[i].velocityGradient[5] = grad_u[1][2];

	  particles[i].velocityGradient[6] = grad_u[2][0];
	  particles[i].velocityGradient[7] = grad_u[2][1];
	  particles[i].velocityGradient[8] = grad_u[2][2];

	  for (j=0;j<prb->nen;j++){
		  particles[i].N0[j] = prb->shape[0][j];
		  particles[i].Nx[j] = prb->shape[1][j*user->iga->dof];
		  particles[i].Ny[j] = prb->shape[1][j*user->iga->dof+1];
		  particles[i].Nz[j] = prb->shape[1][j*user->iga->dof+2];
		  particles[i].map[j] = prb->map[j];
	  }

	  ierr = IGAProbeDestroy(&prb);CHKERRQ(ierr);

  }

  return 0;
}


PetscErrorCode ComputeLumpedMass(void *ctx,PARTICLES *particles,Vec Mass)
{
  PetscErrorCode ierr;
  AppCtx *user = (AppCtx *)ctx;
  PetscInt i,a;
  PetscReal value = 0.0;

  ierr = VecZeroEntries(Mass);CHKERRQ(ierr);

  for (i=0;i<user->numParticles;i++){

	  for (a=0;a<user->nen;a++){
		  PetscInt GlobalID = particles[i].map[a];
		  value = user->Alpha_m*user->rho*particles[i].N0[a]*particles[i].nodalVolume;
		  ierr = VecSetValueLocal(Mass,GlobalID*user->iga->dof,value,ADD_VALUES);CHKERRQ(ierr);
		  ierr = VecSetValueLocal(Mass,GlobalID*user->iga->dof+1,value,ADD_VALUES);CHKERRQ(ierr);
		  ierr = VecSetValueLocal(Mass,GlobalID*user->iga->dof+2,value,ADD_VALUES);CHKERRQ(ierr);
	  }
  }

  ierr = VecAssemblyBegin(Mass);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(Mass);CHKERRQ(ierr);

  return 0;
}



PetscErrorCode GetStress(void *ctx,PARTICLES *particles)
{
  AppCtx *user = (AppCtx *)ctx;
  PetscInt i,j,k,p;

  PetscReal mu = user->mu;
  PetscReal lambda = user->lambda;

  PetscScalar Stress[9] = {0.0};
  PetscScalar StressDev[9] = {0.0};
  PetscScalar DrateDef[9] = {0.0};
  PetscScalar DeltaEpsilon[9] = {0.0};
  PetscScalar DeltaEpsilonDev[9] = {0.0};
  PetscScalar StressTrial[9] = {0.0};
  PetscScalar StressTrialEq = 0.0;
  PetscScalar DeltaEpsilonPlastEq = 0.0;
  PetscScalar stressY = 0.0;
  PetscScalar W[9] = {0.0}; //Spin tensor
  PetscScalar JaumannRate[9] = {0.0};


  for (p=0;p<user->numParticles;p++){

	  Stress[0] = particles[p].Stress0[0];
	  Stress[1] = particles[p].Stress0[3];
	  Stress[2] = particles[p].Stress0[5];
	  Stress[3] = particles[p].Stress0[3];
	  Stress[4] = particles[p].Stress0[1];
	  Stress[5] = particles[p].Stress0[4];
	  Stress[6] = particles[p].Stress0[5];
	  Stress[7] = particles[p].Stress0[4];
	  Stress[8] = particles[p].Stress0[2];

	  for (i=0;i<3;i++){
		  for (j=0;j<3;j++){
			  DrateDef[i*3+j] = 0.5 * (particles[p].velocityGradient[i*3+j] + particles[p].velocityGradient[j*3+i]);
			  W[i*3+j] = 0.5 * (particles[p].velocityGradient[i*3+j] - particles[p].velocityGradient[j*3+i]);
		  }
	  }

	  for (i=0;i<9;i++) JaumannRate[i] = 0.0;

	  for (i=0;i<3;i++){
		  for (j=0;j<3;j++){
			  for (k=0;k<3;k++){
				  JaumannRate[i*3+j] += W[i*3+k]*Stress[k*3+j] + Stress[i*3+k]*W[j*3+k];
			  }
		  }
	  }

	  for(i=0;i<9;i++){
		  Stress[i] = user->Alpha_f * JaumannRate[i] * user->dt + Stress[i];
	  }


	  for (i=0;i<9;i++){
		  DeltaEpsilon[i] = DrateDef[i] * user->dt;
	  }


	  for(i=0;i<9;i++){
		  StressDev[i] = Stress[i];
		  DeltaEpsilonDev[i] = DeltaEpsilon[i];
	  }

	  for(i=0;i<3;i++){
		  StressDev[i*3+i] -= (Stress[0] + Stress[4] + Stress[8]) / 3.0;
		  DeltaEpsilonDev[i*3+i] -= (DeltaEpsilon[0] + DeltaEpsilon[4] + DeltaEpsilon[8]) / 3.0;
	  }

	  for(i=0;i<9;i++){
		  StressTrial[i] = StressDev[i] + 2.0 * mu * DeltaEpsilonDev[i];
	  }


	  StressTrialEq = 0.0;
	  for(i=0;i<9;i++){
		  StressTrialEq += StressTrial[i] * StressTrial[i];
	  }
	  StressTrialEq *= 3.0/2.0;
	  StressTrialEq  = sqrt(StressTrialEq);

	  PetscReal Pressure = (Stress[0] + Stress[4] + Stress[8]) / 3.0;
	  Pressure += (3.0 * lambda + 2.0 * mu) * (DeltaEpsilon[0] + DeltaEpsilon[4] + DeltaEpsilon[8]) / 3.0;

	  stressY = user->stressY + user->h * particles[p].EpsilonPlastEq;

	  if (StressTrialEq - stressY <= 0){
		  for (i=0;i<9;i++) StressDev[i] = StressTrial[i];
	  }else{
		  DeltaEpsilonPlastEq = (StressTrialEq - stressY) / (3.0 * mu + user->h);
		  for (i=0;i<9;i++) StressDev[i] = StressTrial[i] * (1.0 - 3.0 * mu * DeltaEpsilonPlastEq / StressTrialEq);
		  particles[p].DeltaEpsilonPlastEq = DeltaEpsilonPlastEq;
	  }

	  particles[p].Stress[0] = StressDev[0] + Pressure;
	  particles[p].Stress[1] = StressDev[4] + Pressure;
	  particles[p].Stress[2] = StressDev[8] + Pressure;
	  particles[p].Stress[3] = StressDev[1];
	  particles[p].Stress[4] = StressDev[5];
	  particles[p].Stress[5] = StressDev[2];
  }

  return 0;
}



PetscErrorCode ComputeSolidResidual(void *ctx,PARTICLES *particles,Vec Residual)
{
  PetscErrorCode ierr;
  AppCtx *user = (AppCtx *)ctx;
  PetscInt i,a,j;
  PetscReal value = 0.0;
  PetscReal rho = user->rho;
  PetscScalar v[3];
  PetscScalar V[user->nen*user->iga->dof];
  Vec localV;
  const PetscScalar *arrayV;
  ierr = IGAGetLocalVecArray(user->iga,user->Va,&localV,&arrayV);CHKERRQ(ierr);

  ierr = VecZeroEntries(Residual);CHKERRQ(ierr);

  for (i=0;i<user->numParticles;i++){

	  for (a=0;a<user->nen;a++){
		  for (j=0;j<user->iga->dof;j++){
			  V[a*user->iga->dof + j] = arrayV[particles[i].map[a]*user->iga->dof + j];
		  }
	  }

	  v[0] = 0.0;
	  v[1] = 0.0;
	  v[2] = 0.0;
	  for (a=0;a<user->nen;a++){
		  v[0] += particles[i].N0[a] * V[a*user->iga->dof];
		  v[1] += particles[i].N0[a] * V[a*user->iga->dof+1];
		  v[2] += particles[i].N0[a] * V[a*user->iga->dof+2];
	  }

//	  PetscPrintf(PETSC_COMM_WORLD," %d %d \n",prbU->nen,prbV->nen);
	  for (a=0;a<user->nen;a++){
		  PetscInt GlobalID = particles[i].map[a];

		  value = particles[i].Nx[a] * particles[i].Stress[0];
		  value += particles[i].Ny[a] * particles[i].Stress[3];
		  value += particles[i].Nz[a] * particles[i].Stress[5];
		  value += particles[i].N0[a] * rho * v[0];
		  value += particles[i].N0[a] * rho * particles[i].Velocity[0] * particles[i].velocityGradient[0];
		  value += particles[i].N0[a] * rho * particles[i].Velocity[1] * particles[i].velocityGradient[1];
		  value += particles[i].N0[a] * rho * particles[i].Velocity[2] * particles[i].velocityGradient[2];
		  value *= particles[i].nodalVolume;
		  ierr = VecSetValueLocal(Residual,GlobalID*user->iga->dof,value,ADD_VALUES);CHKERRQ(ierr);

		  value = particles[i].Nx[a] * particles[i].Stress[3];
		  value += particles[i].Ny[a] * particles[i].Stress[1];
		  value += particles[i].Nz[a] * particles[i].Stress[4];
		  value += particles[i].N0[a] * rho * v[1];
		  value += particles[i].N0[a] * rho * particles[i].Velocity[0] * particles[i].velocityGradient[3];
		  value += particles[i].N0[a] * rho * particles[i].Velocity[1] * particles[i].velocityGradient[4];
		  value += particles[i].N0[a] * rho * particles[i].Velocity[2] * particles[i].velocityGradient[5];
		  value *= particles[i].nodalVolume;
		  ierr = VecSetValueLocal(Residual,GlobalID*user->iga->dof+1,value,ADD_VALUES);CHKERRQ(ierr);

		  value = particles[i].Nx[a] * particles[i].Stress[5];
		  value += particles[i].Ny[a] * particles[i].Stress[4];
		  value += particles[i].Nz[a] * particles[i].Stress[2];
		  value += particles[i].N0[a] * rho * v[2];
		  value += particles[i].N0[a] * rho * particles[i].Velocity[0] * particles[i].velocityGradient[6];
		  value += particles[i].N0[a] * rho * particles[i].Velocity[1] * particles[i].velocityGradient[7];
		  value += particles[i].N0[a] * rho * particles[i].Velocity[2] * particles[i].velocityGradient[8];
		  value *= particles[i].nodalVolume;
		  ierr = VecSetValueLocal(Residual,GlobalID*user->iga->dof+2,value,ADD_VALUES);CHKERRQ(ierr);
	  }


  }

  ierr = VecAssemblyBegin(Residual);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(Residual);CHKERRQ(ierr);

  ierr = IGARestoreLocalVecArray(user->iga,user->Va,&localV,&arrayV);CHKERRQ(ierr);

  return 0;
}


PetscErrorCode SetSmallMass(void *ctx,Vec Mass,Vec F)
{
  PetscErrorCode ierr;
  AppCtx *user = (AppCtx *)ctx;
  PetscInt m,l,k;

  PetscInt dof = user->iga->dof;

  PetscScalar value[1];
  PetscInt index2[1];
  PetscInt index;

  PetscInt nodesX  = user->iga->geom_lwidth[0];
  PetscInt nodesY  = user->iga->geom_lwidth[1];
  PetscInt nodesZ  = user->iga->geom_lwidth[2];

  for (m=0;m<nodesZ;m++){
      for (l=0;l<nodesY;l++){
          for (k=0;k<nodesX;k++){
              index2[0] = (m*nodesX*nodesY + l*nodesX + k)*dof;
              ierr = VecGetValues(Mass,1,index2,value);CHKERRQ(ierr);

              ierr = VecAssemblyBegin(Mass);CHKERRQ(ierr);
              ierr = VecAssemblyEnd(Mass);CHKERRQ(ierr);

              if (value[0]<1.0e-10){
            	  index = (m*nodesX*nodesY + l*nodesX + k)*dof;
            	  VecSetValue(Mass,index,1.0,INSERT_VALUES);CHKERRQ(ierr);
    	          VecSetValue(F,index,0.0,INSERT_VALUES);CHKERRQ(ierr);

    	          index = (m*nodesX*nodesY + l*nodesX + k)*dof+1;
    	          VecSetValue(Mass,index,1.0,INSERT_VALUES);CHKERRQ(ierr);
    	          VecSetValue(F,index,0.0,INSERT_VALUES);CHKERRQ(ierr);

    	          index = (m*nodesX*nodesY + l*nodesX + k)*dof+2;
    	          VecSetValue(Mass,index,1.0,INSERT_VALUES);CHKERRQ(ierr);
    	          VecSetValue(F,index,0.0,INSERT_VALUES);CHKERRQ(ierr);
              }

          }
      }
  }

  ierr = VecAssemblyBegin(Mass);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(Mass);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F);CHKERRQ(ierr);

  return 0;
}



PetscErrorCode SetBoundaryConditions(Vec Mass,Vec F,void *ctx)
{
  PetscErrorCode ierr;
  AppCtx *user = (AppCtx *)ctx;

  PetscInt m,l,k;
  PetscInt dim = user->iga->dim;
  PetscInt dof = user->iga->dof;

  PetscInt nodesX  = user->iga->geom_lwidth[0];
  PetscInt nodesY  = user->iga->geom_lwidth[1];
  PetscInt nodesZ  = user->iga->geom_lwidth[2];
  PetscInt gnodesX = user->iga->geom_gwidth[0];
  PetscInt gnodesY = user->iga->geom_gwidth[1];

  for (m=0;m<nodesZ;m++){
      for (l=0;l<nodesY;l++){
          for (k=0;k<nodesX;k++){

	            if (user->iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX + k)*dim + 2] <= 0.00001){
	            	PetscInt index = (m*gnodesX*gnodesY + l*gnodesX + k)*dof+2;
	            	ierr = VecSetValueLocal(Mass,index,1.0,INSERT_VALUES);CHKERRQ(ierr);
	            	ierr = VecSetValueLocal(F,index,0.0,INSERT_VALUES);CHKERRQ(ierr);
	            }

//	            if (user->iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX + k)*dim] <= 0.00001){
//	            	PetscInt index = (m*gnodesX*gnodesY + l*gnodesX + k)*dof;
//	            	ierr = VecSetValueLocal(Mass,index,1.0,INSERT_VALUES);CHKERRQ(ierr);
//	            	ierr = VecSetValueLocal(F,index,0.0,INSERT_VALUES);CHKERRQ(ierr);
//	            }
//
//	            if (user->iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX + k)*dim + 1] <= 0.00001){
//	            	PetscInt index = (m*gnodesX*gnodesY + l*gnodesX + k)*dof+1;
//	            	ierr = VecSetValueLocal(Mass,index,1.0,INSERT_VALUES);CHKERRQ(ierr);
//	            	ierr = VecSetValueLocal(F,index,0.0,INSERT_VALUES);CHKERRQ(ierr);
//	            }

          }
      }
  }

  ierr = VecAssemblyBegin(Mass);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(Mass);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F);CHKERRQ(ierr);

  return 0;
}


PetscErrorCode SolveSystem(Vec Mass,Vec F,void *ctx)
{
  PetscErrorCode ierr;
  AppCtx *user = (AppCtx *)ctx;

  ierr = VecReciprocal(Mass);CHKERRQ(ierr);
  ierr = VecPointwiseMult(user->dV,Mass,F);CHKERRQ(ierr);

  return 0;
}


PetscErrorCode UpdateStage(void *ctx)
{
  PetscErrorCode ierr;
  AppCtx *user = (AppCtx *)ctx;

  ierr = VecAXPY(user->V1,-1.0,user->dV);CHKERRQ(ierr);
  ierr = VecAXPY(user->U1,-user->Gamma*user->dt,user->dV);CHKERRQ(ierr);

  return 0;
}


PetscErrorCode InterpolateVelocityOnParticles(void *ctx,PARTICLES *particles)
{
  PetscErrorCode ierr;
  AppCtx *user = (AppCtx *)ctx;
  PetscInt i;

  PetscScalar u[3] = {0.0};
  PetscScalar pt[3] = {0.0};


  for (i=0;i<user->numParticles;i++){

	  pt[0] = particles[i].currentCoord[0] / user->Lx;
	  pt[1] = particles[i].currentCoord[1] / user->Ly;
	  pt[2] = particles[i].currentCoord[2] / user->Lz;

	  IGAProbe prb;
	  ierr = IGAProbeCreate(user->iga,user->U1,&prb);CHKERRQ(ierr);
	  ierr = IGAProbeSetPoint(prb,pt);CHKERRQ(ierr);
	  ierr = IGAProbeFormValue(prb,&u[0]);CHKERRQ(ierr);
	  ierr = IGAProbeDestroy(&prb);CHKERRQ(ierr);

	  particles[i].Velocity[0] = u[0];
	  particles[i].Velocity[1] = u[1];
	  particles[i].Velocity[2] = u[2];

  }

  return 0;
}


PetscErrorCode UpdateParticles(void *ctx,PARTICLES *particles)
{
  PetscErrorCode ierr;
  AppCtx *user = (AppCtx *)ctx;
  PetscInt i;

  ierr = InterpolateVelocityOnParticles(user,particles);CHKERRQ(ierr);

  for (i=0;i<user->numParticles;i++){
	  particles[i].AccelerationIncrement[0] = particles[i].Velocity[0] - particles[i].VelocityOldIteration[0];
	  particles[i].AccelerationIncrement[0] /= user->Gamma * user->dt;

	  particles[i].AccelerationIncrement[1] = particles[i].Velocity[1] - particles[i].VelocityOldIteration[1];
	  particles[i].AccelerationIncrement[1] /= user->Gamma * user->dt;

	  particles[i].AccelerationIncrement[2] = particles[i].Velocity[2] - particles[i].VelocityOldIteration[2];
	  particles[i].AccelerationIncrement[2] /= user->Gamma * user->dt;


	  particles[i].Acceleration[0] = particles[i].AccelerationOldIteration[0] + particles[i].AccelerationIncrement[0];
	  particles[i].Acceleration[1] = particles[i].AccelerationOldIteration[1] + particles[i].AccelerationIncrement[1];
	  particles[i].Acceleration[2] = particles[i].AccelerationOldIteration[2] + particles[i].AccelerationIncrement[2];

	  particles[i].Displacement[0] = particles[i].DisplacementOldIteration[0];
	  particles[i].Displacement[0] += user->Beta * user->dt * user->dt * particles[i].AccelerationIncrement[0];

	  particles[i].Displacement[1] = particles[i].DisplacementOldIteration[1];
	  particles[i].Displacement[1] += user->Beta * user->dt * user->dt * particles[i].AccelerationIncrement[1];

	  particles[i].Displacement[2] = particles[i].DisplacementOldIteration[2];
	  particles[i].Displacement[2] += user->Beta * user->dt * user->dt * particles[i].AccelerationIncrement[2];

	  particles[i].currentCoord[0] = particles[i].initialCoord[0] + particles[i].Displacement[0];
	  particles[i].currentCoord[1] = particles[i].initialCoord[1] + particles[i].Displacement[1];
	  particles[i].currentCoord[2] = particles[i].initialCoord[2] + particles[i].Displacement[2];
  }

  return 0;
}


PetscErrorCode RotateStress(void *ctx,PARTICLES *particles)
{
  AppCtx *user = (AppCtx *)ctx;
  PetscInt i,j,k,p;

  PetscScalar Stress[9] = {0.0};
  PetscScalar W[9] = {0.0}; //Spin tensor
  PetscScalar JaumannRate[9] = {0.0};


  for (p=0;p<user->numParticles;p++){

	  Stress[0] = particles[p].Stress[0];
	  Stress[1] = particles[p].Stress[3];
	  Stress[2] = particles[p].Stress[5];
	  Stress[3] = particles[p].Stress[3];
	  Stress[4] = particles[p].Stress[1];
	  Stress[5] = particles[p].Stress[4];
	  Stress[6] = particles[p].Stress[5];
	  Stress[7] = particles[p].Stress[4];
	  Stress[8] = particles[p].Stress[2];

	  for (i=0;i<3;i++){
		  for (j=0;j<3;j++){
			  W[i*3+j] = 0.5 * (particles[p].velocityGradient[i*3+j] - particles[p].velocityGradient[j*3+i]);
		  }
	  }

	  for (i=0;i<9;i++) JaumannRate[i] = 0.0;

	  for (i=0;i<3;i++){
		  for (j=0;j<3;j++){
			  for (k=0;k<3;k++){
				  JaumannRate[i*3+j] += W[i*3+k]*Stress[k*3+j] + Stress[i*3+k]*W[j*3+k];
			  }
		  }
	  }

	  for(i=0;i<9;i++){
		  Stress[i] = (1.0 - user->Alpha_f) * JaumannRate[i] * user->dt + Stress[i];
	  }

	  particles[p].Stress0[0] = Stress[0];
	  particles[p].Stress0[1] = Stress[4];
	  particles[p].Stress0[2] = Stress[8];
	  particles[p].Stress0[3] = Stress[1];
	  particles[p].Stress0[4] = Stress[5];
	  particles[p].Stress0[5] = Stress[2];

	  particles[p].EpsilonPlastEq += particles[p].DeltaEpsilonPlastEq;
  }

  return 0;
}


PetscErrorCode OutputInitial(void *ctx,PARTICLES *particles)
{
  PetscErrorCode ierr;
  AppCtx *user = (AppCtx *)ctx;
  FILE *fileResult;
  char filenameResults[50];

  ierr = IGAWrite(user->iga,"Geo.dat");CHKERRQ(ierr);
  ierr = IGAWriteVec(user->iga,user->U0,"Vel0.dat");CHKERRQ(ierr);


  sprintf(filenameResults,"./Meshfree.case");
  fileResult=fopen(filenameResults,"wt");

  PetscInt numSteps = (int) (user->final_time/user->dt);

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
  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"vector per node: 1 Velocity Velocity.*.res \n");
  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"\n");
  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"TIME \n");
  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"\n");

  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"time set: 1 \n");
  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"number of steps: %d \n", numSteps/user->freqResults);
  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"filename start number: 0 \n");
  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"filename increment: %d \n", 1);
  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"time values: \n");

  PetscInt counter = 0;
  PetscInt i;
  for (i=0;i<numSteps;i++){
	  if (i % user->freqResults == 0){
		  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"%d ",counter);
		  counter++;
	  }
	  if ((i+1) % (10*user->freqResults) == 0) PetscFPrintf(PETSC_COMM_WORLD,fileResult,"\n  ");
  }

  fclose(fileResult);



  sprintf(filenameResults,"./Meshless.0.geo");
  fileResult=fopen(filenameResults,"wt");
  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"Meshless \n");
  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"node \n");
  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"node id given \n");
  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"element id given \n");
  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"coordinates \n");
  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"%d \n",user->numParticles);

  for (i=0;i<user->numParticles;i++){
	  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"%d   %.4e   %.4e   %.4e \n",i+1,
			  particles[i].currentCoord[0],
			  particles[i].currentCoord[1],
			  particles[i].currentCoord[2]);

  }

  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"part     1\n");
  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"todo \n");
  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"point \n");

  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"%d \n",user->numParticles);

  for(i=1;i<=user->numParticles;i++){
  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"%d  %d \n",i,i);
  }

  fclose(fileResult);



  sprintf(filenameResults,"./Velocity.0.res");
  fileResult=fopen(filenameResults,"wt");
  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"Velocity \n");

  for (i=0;i<user->numParticles;i++){
	  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"%.4e   %.4e   %.4e  ",particles[i].Velocity[0],
			  particles[i].Velocity[1],
			  particles[i].Velocity[2]);

	  if ((i+1) % 2 == 0) PetscFPrintf(PETSC_COMM_WORLD,fileResult,"\n  ");
  }

  fclose(fileResult);


  return 0;
}


PetscErrorCode Output(void *ctx,PARTICLES *particles)
{
  PetscErrorCode ierr;
  AppCtx *user = (AppCtx *)ctx;
  PetscInt step = user->step;
  PetscInt i;
  PetscPrintf(PETSC_COMM_WORLD," Step %d \n",step);
  if (step % user->freqResults == 0){
	  char filename[256];
	  sprintf(filename,"Vel%d.dat",step/user->freqResults);
	  ierr = IGAWriteVec(user->iga,user->U1,filename);CHKERRQ(ierr);


	  FILE *fileResult;
	   char filenameResults[50];

	   sprintf(filenameResults,"./Meshless.%d.geo",step/user->freqResults);
	   fileResult=fopen(filenameResults,"wt");
	   PetscFPrintf(PETSC_COMM_WORLD,fileResult,"Meshless \n");
	   PetscFPrintf(PETSC_COMM_WORLD,fileResult,"node \n");
	   PetscFPrintf(PETSC_COMM_WORLD,fileResult,"node id given \n");
	   PetscFPrintf(PETSC_COMM_WORLD,fileResult,"element id given \n");
	   PetscFPrintf(PETSC_COMM_WORLD,fileResult,"coordinates \n");
	   PetscFPrintf(PETSC_COMM_WORLD,fileResult,"%d \n",user->numParticles);

	   for (i=0;i<user->numParticles;i++){
	 	  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"%d   %.4e   %.4e   %.4e \n",i+1,
	 			  particles[i].currentCoord[0],
	 			  particles[i].currentCoord[1],
	 			  particles[i].currentCoord[2]);

	   }

	   PetscFPrintf(PETSC_COMM_WORLD,fileResult,"part     1\n");
	   PetscFPrintf(PETSC_COMM_WORLD,fileResult,"todo \n");
	   PetscFPrintf(PETSC_COMM_WORLD,fileResult,"point \n");

	   PetscFPrintf(PETSC_COMM_WORLD,fileResult,"%d \n",user->numParticles);

	   for(i=1;i<=user->numParticles;i++){
	   PetscFPrintf(PETSC_COMM_WORLD,fileResult,"%d  %d \n",i,i);
	   }

	   fclose(fileResult);



	   sprintf(filenameResults,"./Velocity.%d.res",step/user->freqResults);
	   fileResult=fopen(filenameResults,"wt");
	   PetscFPrintf(PETSC_COMM_WORLD,fileResult,"Velocity \n");

	   for (i=0;i<user->numParticles;i++){
	 	  PetscFPrintf(PETSC_COMM_WORLD,fileResult,"%.4e   %.4e   %.4e  ",particles[i].Velocity[0],
	 			  particles[i].Velocity[1],
	 			  particles[i].Velocity[2]);

	 	  if ((i+1) % 2 == 0) PetscFPrintf(PETSC_COMM_WORLD,fileResult,"\n  ");
	   }

	   fclose(fileResult);


  }

  return 0;
}


PetscErrorCode UpdateStep(void *ctx)
{
  PetscErrorCode ierr;
  AppCtx *user = (AppCtx *)ctx;

  ierr = VecCopy(user->U1,user->U0);CHKERRQ(ierr);
  ierr = VecCopy(user->V1,user->V0);CHKERRQ(ierr);

  return 0;
}


PetscErrorCode SetRadius_GeneralizedAlpha(AppCtx *user,PetscReal radius)
{
  user->Alpha_m = (3.0-radius)/(2.0*(1.0 + radius));
  user->Alpha_f = 1.0/(1.0+radius);
  user->Gamma   = 0.5 + user->Alpha_m - user->Alpha_f;
  user->Beta    = 0.5 * (1.0 + user->Alpha_m - user->Alpha_f);
  user->Beta *= user->Beta;
  return 0;
}

int main(int argc, char *argv[]) {

  PetscErrorCode ierr;
  ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);


  AppCtx user;
  user.E = 78.2E+09; //Young's modulus
  user.nu = 0.3; //Poisson ratio
  user.rho = 2700.0; //density
  user.lambda = user.E * user.nu / ((1.0 + user.nu) * (1.0 - 2.0 * user.nu)); //Lame parameter
  user.mu = user.E / (2 * (1 + user.nu)); //Lame parameter
  user.dt = 1.0e-8; //time step size
  user.final_time = 4.0e-5; //analysis time
  user.stressY = 0.29e+09; //initial yield stress
  user.h = 1.0e+8; //hardening parameter
  user.maxIt = 3; //number of iterations
  user.freqResults = 100;


  PetscInt dim; //spatial dimensions
  IGA iga;
  ierr = IGACreate(PETSC_COMM_WORLD,&iga);CHKERRQ(ierr);
  ierr = IGARead(iga,"TaylorGeometry.dat");CHKERRQ(ierr);
  ierr = IGAGetDim(iga,&dim);CHKERRQ(ierr);
  ierr = IGASetDof(iga,dim);CHKERRQ(ierr);
  ierr = IGASetOrder(iga,2);CHKERRQ(ierr);
  ierr = IGASetFromOptions(iga);CHKERRQ(ierr);
  ierr = IGASetUp(iga);CHKERRQ(ierr);


  //U velocity
  //V acceleration
  //0 - zero level (n level)
  //p - predictor
  //a - alpha level
  //1 - n+1 level

  ierr = IGACreateVec(iga,&user.dV);CHKERRQ(ierr);
  ierr = VecZeroEntries(user.dV);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&user.U0);CHKERRQ(ierr);
  ierr = VecZeroEntries(user.U0);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&user.Up);CHKERRQ(ierr);
  ierr = VecZeroEntries(user.Up);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&user.U1);CHKERRQ(ierr);
  ierr = VecZeroEntries(user.U1);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&user.Ua);CHKERRQ(ierr);
  ierr = VecZeroEntries(user.Ua);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&user.V0);CHKERRQ(ierr);
  ierr = VecZeroEntries(user.V0);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&user.Vp);CHKERRQ(ierr);
  ierr = VecZeroEntries(user.Vp);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&user.V1);CHKERRQ(ierr);
  ierr = VecZeroEntries(user.V1);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&user.Va);CHKERRQ(ierr);
  ierr = VecZeroEntries(user.Va);CHKERRQ(ierr);

  Vec F; //Force Vector
  Vec Mass; //Lumped mass vector
  ierr = IGACreateVec(iga,&F);CHKERRQ(ierr);
  ierr = VecZeroEntries(F);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&Mass);CHKERRQ(ierr);
  ierr = VecZeroEntries(Mass);CHKERRQ(ierr);


  user.iga = iga;
  ierr = SetRadius_GeneralizedAlpha(&user,0.5);CHKERRQ(ierr);

  PARTICLES *particles;
  ierr = ReadNumParticles(&user);CHKERRQ(ierr);
  ierr = PetscMalloc1(sizeof(*particles)*user.numParticles,&particles);CHKERRQ(ierr);
  ierr = Input(&user,particles);CHKERRQ(ierr);

  PetscReal t = 0.0; //current time
  user.step = 0; //current step
  ierr = FormInitialCondition(&user,iga,t,user.U0);CHKERRQ(ierr);
  ierr = OutputInitial(&user,particles);CHKERRQ(ierr);

//  VecView(U,PETSC_VIEWER_STDOUT_WORLD);

  while (t < user.final_time){

	  PetscInt i;
	  t += user.dt;
	  user.step += 1;
	  user.it = 0;

	  ierr = SetParticlesOldStep(&user,particles);CHKERRQ(ierr);
	  ierr = PredictorStage(&user,particles);CHKERRQ(ierr);

	  for (i=0;i<user.maxIt;i++){

		  user.it += 1;
		  ierr = SetParticlesOldIteration(&user,particles);CHKERRQ(ierr);
		  ierr = AlphaStage(&user,particles);CHKERRQ(ierr);
		  ierr = ComputeVelocityGradientAndShapeFunctions(&user,particles);CHKERRQ(ierr);
		  ierr = ComputeLumpedMass(&user,particles,Mass);CHKERRQ(ierr);
		  ierr = GetStress(&user,particles);CHKERRQ(ierr);
		  ierr = ComputeSolidResidual(&user,particles,F);CHKERRQ(ierr);
		  ierr = SetSmallMass(&user,Mass,F);CHKERRQ(ierr);
		  ierr = SetBoundaryConditions(Mass,F,&user);CHKERRQ(ierr);
		  ierr = SolveSystem(Mass,F,&user);CHKERRQ(ierr);
		  ierr = UpdateStage(&user);CHKERRQ(ierr);
		  ierr = UpdateParticles(&user,particles);CHKERRQ(ierr);

	  } //End of iteration loop

//	  ierr = ComputeVelocityGradient(&user,particles);CHKERRQ(ierr);
	  ierr = RotateStress(&user,particles);CHKERRQ(ierr);
	  ierr = Output(&user,particles);CHKERRQ(ierr);
	  ierr = UpdateStep(&user);CHKERRQ(ierr);

  } //End of time loop


  ierr = IGADestroy(&iga);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}
