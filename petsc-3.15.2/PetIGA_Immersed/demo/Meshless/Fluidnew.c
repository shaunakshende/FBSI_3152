#include "petiga.h"
#include "time.h"


typedef struct {
 PetscReal lamda,Lx,Ly,Lz,kappa,Cp,R,Cv,mu;
 PetscScalar F[4];
 IGA iga;
 PetscReal temp0;
 PetscReal p0;
 PetscReal currentTime;
 PetscReal timeStep;
 PetscReal finalTime;
 PetscInt  stepNumber;

 PetscReal Alpha_m,Alpha_f,Gamma,Beta;


 Vec V0,A0;
 Vec V1,A1;
 Vec Vp,Ap;
 Vec Va,Aa;
 Vec dA;

 PetscInt max_its;

 PetscReal TimeRestart;
 PetscInt  StepRestart;
 PetscInt  FreqRestarts;
 PetscInt  FreqResults;

} AppCtx;

#undef __FUNCT__
#define __FUNCT__ "Tau"
PetscErrorCode Tau(PetscReal G[2][2],
         PetscReal dt,PetscScalar u[],
         PetscScalar (*tau)[4],void *ctx,PetscReal (*A0inv)[4],PetscReal *umi)
{

	AppCtx *user = (AppCtx *)ctx;

	 PetscInt  dim    = user->iga->dim;
	 PetscInt  dof    = user->iga->dof;
     ;
     PetscReal kappa  = user->kappa;
     PetscReal Cv     = user->Cv;
     PetscReal RR     = user->R;
     PetscReal P0 = 100000;
     PetscReal B  = 3000*P0;
     PetscReal N  = 7.14;

     PetscInt i,j,l;
	   PetscScalar P    = u[0];
     PetscScalar ux   = u[1];
     PetscScalar uy   = u[2];
     PetscScalar temp = u[3];
     PetscReal dens0    = 1000;
     PetscReal dens     = dens0*pow((1/B)*(P-P0)+1, 1/N);
     PetscReal f        = P0+B*(pow(dens/dens0,N)-1);
     PetscReal fprime   = (1.0/dens)*B*N*(pow(dens/dens0,N));
     PetscReal fcr        = P0+B*(pow(0.99,N)-1);
     if(P<=fcr+0.001){
       dens = dens0*pow((1/B)*(fcr-P0)+1, 1/N);
     }
     f        = P0+B*(pow(dens/dens0,N)-1);
     fprime   = (1.0/dens)*B*N*(pow(dens/dens0,N));

     PetscReal mu       = user->mu;

     PetscReal Ad_vel_l[2];
     Ad_vel_l[0]=ux-umi[0];
     Ad_vel_l[1]=uy-umi[1];

     PetscReal gamma_c=1.4;
     PetscReal cs=sqrt(fprime);
//
////PetscPrintf(PETSC_COMM_WORLD,"P %e \n",P );
////PetscPrintf(PETSC_COMM_WORLD,"ux %e \n",ux );
////PetscPrintf(PETSC_COMM_WORLD,"uy %e \n",uy );
////PetscPrintf(PETSC_COMM_WORLD,"uz %e \n",uz );
////PetscPrintf(PETSC_COMM_WORLD,"temp %e \n",temp );
////PetscPrintf(PETSC_COMM_WORLD,"dens %e \n",dens );
////PetscPrintf(PETSC_COMM_WORLD,"P %e \n",P );
////PetscPrintf(PETSC_COMM_WORLD,"RR %e \n",RR );
//
     PetscReal Evec1[2]={0.0};
     PetscReal Evec2[2]={0.0};
     Evec1[0]=1.0;
     Evec2[1]=1.0;

     PetscReal Umag=sqrt(Ad_vel_l[0]*Ad_vel_l[0]+Ad_vel_l[1]*Ad_vel_l[1]);
	//PetscPrintf(PETSC_COMM_WORLD,"Umag %e \n",Umag );
   if (Umag>=1e-9){
	   Evec1[0]=Ad_vel_l[0]/Umag;
	   Evec1[1]=Ad_vel_l[1]/Umag;

	   Evec2[0]=-Ad_vel_l[1]/Umag;
	   Evec2[1]= Ad_vel_l[0]/Umag;

////   PetscPrintf(PETSC_COMM_WORLD,"Emag %e \n",Emag );
   }
//
//
//
//
////   PetscPrintf(PETSC_COMM_WORLD,"Evec3[0] %e \n",Evec3[0] );
////   PetscPrintf(PETSC_COMM_WORLD,"Evec3[1] %e \n",Evec3[1] );
////   PetscPrintf(PETSC_COMM_WORLD,"Evec3[2] %e \n",Evec3[2] );
////   PetscPrintf(PETSC_COMM_WORLD,"Evec1[0] %e \n",Evec1[0] );
////   PetscPrintf(PETSC_COMM_WORLD,"Evec1[1] %e \n",Evec1[1] );
////   PetscPrintf(PETSC_COMM_WORLD,"Evec1[2] %e \n",Evec1[2] );
////   PetscPrintf(PETSC_COMM_WORLD,"Evec2[0] %e \n",Evec2[0] );
////   PetscPrintf(PETSC_COMM_WORLD,"Evec2[1] %e \n",Evec2[1] );
////   PetscPrintf(PETSC_COMM_WORLD,"Evec2[2] %e \n",Evec2[2] );
//
   PetscReal alpha1=0.0;
   PetscReal alpha2=0.0;
   for (i=0;i<dim;i++)
	   for (j=0;j<dim;j++){
		   alpha1=alpha1+Evec1[i]*G[i][j]*Evec1[j];
   	   	   alpha2=alpha2+Evec2[i]*G[i][j]*Evec2[j];
	   }



   PetscReal alpha3=0.0;

 // Tau for Conservation variables

   PetscReal taua=0.0;
   PetscReal taut=0.0;
   PetscReal taud=0.0;
   PetscReal gij2=0.0;

   taut = 4.0/(dt*dt);
   taua = alpha1*(Umag*Umag + cs*cs) + 0.5*(alpha2+alpha3)*cs*cs + 0.5*cs*sqrt(cs*cs*(alpha2 + alpha3)*(alpha2 + alpha3) + 16.0*alpha1*alpha1*Umag*Umag);

////   PetscPrintf(PETSC_COMM_WORLD,"alpha1 %e \n",alpha1 );
////   PetscPrintf(PETSC_COMM_WORLD,"alpha2 %e \n",alpha2 );
////   PetscPrintf(PETSC_COMM_WORLD,"alpha3 %e \n",alpha3 );
////   PetscPrintf(PETSC_COMM_WORLD,"Umag %e \n",Umag );
////   PetscPrintf(PETSC_COMM_WORLD,"cs %e \n",cs );
////   PetscPrintf(PETSC_COMM_WORLD,"taua %e \n",taua );
//
  PetscReal tauC=1.0/sqrt(taut+taua);

  for (j=0;j<dim;j++)
	   for (i=0;i<dim;i++)
		   gij2=gij2+G[i][j]*G[i][j];

  PetscReal m_k=3.0;

  taud=m_k*mu*mu*gij2/(dens*dens);
  PetscReal tauM=1.0/sqrt(taut+taua+taud);

  taud=m_k*kappa*kappa*gij2/(dens*dens*Cv*Cv);
  PetscReal tauE=1.0/sqrt(taut+taua+taud);

////  PetscPrintf(PETSC_COMM_WORLD,"taut %e \n",taut );
////  PetscPrintf(PETSC_COMM_WORLD,"taua %e \n",taua );
////  PetscPrintf(PETSC_COMM_WORLD,"tauM %e \n",tauM );
////  PetscPrintf(PETSC_COMM_WORLD,"tauE %e \n",tauE );cd
//
 // Transform to pressure primitive variables

  PetscReal TauTemp[4][4]={{0.0}};
  TauTemp[0][0] = tauC;
  TauTemp[1][1] = tauM;
  TauTemp[2][2] = tauM;
  TauTemp[3][3] = tauE;

  PetscInt aa,bb;
  for (aa=0;aa<dof;aa++)
 	   for (bb=0;bb<dof;bb++)
 		   tau[aa][bb] = A0inv[aa][0]*TauTemp[0][bb]
 		                +A0inv[aa][1]*TauTemp[1][bb]
 		                +A0inv[aa][2]*TauTemp[2][bb]
 		                +A0inv[aa][3]*TauTemp[3][bb];





  return 0;
}





#undef __FUNCT__
#define __FUNCT__ "Tau1"
PetscErrorCode Tau1(PetscReal G[2][2],
         PetscReal dt,PetscScalar u[],
         PetscScalar (*tau)[4],void *ctx,PetscReal (*A0inv)[4],PetscReal *umi)
{

	AppCtx *user = (AppCtx *)ctx;

	 PetscInt  dim    = user->iga->dim;
	 PetscInt  dof    = user->iga->dof;
     PetscReal kappa  = user->kappa;
     PetscReal Cv     = user->Cv;
     //PetscReal mu     = user->mu;
     PetscReal RR     = user->R;
     PetscReal B = 3000*100000;
     PetscReal P0 = 100000;
     PetscReal N = 7.14;

     PetscInt i,j;
	   PetscScalar P    = u[0];
     PetscScalar ux   = u[1];
     PetscScalar uy   = u[2];
     PetscScalar temp = u[3];
     PetscReal dens0    = 1000;
     PetscReal dens     = dens0*pow((1/B)*(P-P0)+1, 1/N);
     PetscReal f        = P0+B*(pow(dens/dens0,N)-1);
     PetscReal fprime   = (1.0/dens)*B*N*(pow(dens/dens0,N));
     PetscReal fcr        = P0+B*(pow(0.99,N)-1);
     if(P<=fcr+0.001){
       dens = dens0*pow((1/B)*(fcr-P0)+1, 1/N);
     }
     f        = P0+B*(pow(dens/dens0,N)-1);
     fprime   = (1.0/dens)*B*N*(pow(dens/dens0,N));

     PetscReal mu       = user->mu;

     PetscReal Ad_vel_l[2];
     Ad_vel_l[0]=ux-umi[0];
     Ad_vel_l[1]=uy-umi[1];

     PetscReal gamma_c=1.4;
     PetscReal cs=sqrt(fprime);

//PetscPrintf(PETSC_COMM_WORLD,"P %e \n",P );
//PetscPrintf(PETSC_COMM_WORLD,"ux %e \n",ux );
//PetscPrintf(PETSC_COMM_WORLD,"uy %e \n",uy );
//PetscPrintf(PETSC_COMM_WORLD,"uz %e \n",uz );
//PetscPrintf(PETSC_COMM_WORLD,"temp %e \n",temp );
//PetscPrintf(PETSC_COMM_WORLD,"dens %e \n",dens );
//PetscPrintf(PETSC_COMM_WORLD,"P %e \n",P );
//PetscPrintf(PETSC_COMM_WORLD,"RR %e \n",RR );

     PetscReal Evec1[2]={0.0};
     PetscReal Evec2[2]={0.0};
     Evec1[0]=1.0;
     Evec2[1]=1.0;

     PetscReal Umag=sqrt(Ad_vel_l[0]*Ad_vel_l[0]+Ad_vel_l[1]*Ad_vel_l[1]);


     PetscReal tauC_bar = 0.5*dt;
     PetscReal tauM_bar = 0.5*dt;
     PetscReal tauE_bar = 0.5*dt;

   if (Umag>=1e-9){
	   Evec1[0]=Ad_vel_l[0]/Umag;
	   Evec1[1]=Ad_vel_l[1]/Umag;

	   Evec2[0]=-Ad_vel_l[1]/Umag;
	   Evec2[1]= Ad_vel_l[0]/Umag;

//   PetscPrintf(PETSC_COMM_WORLD,"Emag %e \n",Emag );




//   PetscPrintf(PETSC_COMM_WORLD,"Evec3[0] %e \n",Evec3[0] );
//   PetscPrintf(PETSC_COMM_WORLD,"Evec3[1] %e \n",Evec3[1] );
//   PetscPrintf(PETSC_COMM_WORLD,"Evec3[2] %e \n",Evec3[2] );
//   PetscPrintf(PETSC_COMM_WORLD,"Evec1[0] %e \n",Evec1[0] );
//   PetscPrintf(PETSC_COMM_WORLD,"Evec1[1] %e \n",Evec1[1] );
//   PetscPrintf(PETSC_COMM_WORLD,"Evec1[2] %e \n",Evec1[2] );
//   PetscPrintf(PETSC_COMM_WORLD,"Evec2[0] %e \n",Evec2[0] );
//   PetscPrintf(PETSC_COMM_WORLD,"Evec2[1] %e \n",Evec2[1] );
//   PetscPrintf(PETSC_COMM_WORLD,"Evec2[2] %e \n",Evec2[2] );

   PetscReal alpha1=0.0;

   for (i=0;i<dim;i++)
	   for (j=0;j<dim;j++)
		   alpha1=alpha1+Evec1[i]*G[i][j]*Evec1[j];


   PetscReal alpha2=0.0;

   for (i=0;i<dim;i++)
	   for (j=0;j<dim;j++)
		   alpha2=alpha2+Evec2[i]*G[i][j]*Evec2[j];

   PetscReal alpha3=0.0;

  //Tau for Conservation variables

   PetscReal taua=0.0;
   PetscReal taut=0.0;
   PetscReal taud=0.0;
   PetscReal gij2=0.0;

   taut = 4.0/(dt*dt);
   taua = alpha1*(Umag*Umag + cs*cs) + 0.5*(alpha2+alpha3)*cs*cs + 0.5*cs*sqrt(cs*cs*(alpha2 + alpha3)*(alpha2 + alpha3) + 16.0*alpha1*alpha1*Umag*Umag);

//   PetscPrintf(PETSC_COMM_WORLD,"alpha1 %e \n",alpha1 );
//   PetscPrintf(PETSC_COMM_WORLD,"alpha2 %e \n",alpha2 );
//   PetscPrintf(PETSC_COMM_WORLD,"alpha3 %e \n",alpha3 );
//   PetscPrintf(PETSC_COMM_WORLD,"Umag %e \n",Umag );
//   PetscPrintf(PETSC_COMM_WORLD,"cs %e \n",cs );
//   PetscPrintf(PETSC_COMM_WORLD,"taua %e \n",taua );

  PetscReal tauC=1.0/sqrt(taut+taua);



  for (j=0;j<dim;j++)
	   for (i=0;i<dim;i++)
		   gij2=gij2+G[i][j]*G[i][j];

  PetscReal m_k=3.0;

  taud=m_k*mu*mu*gij2/(dens*dens);
  PetscReal tauM=1.0/sqrt(taut+taua+taud);

  taud=m_k*kappa*kappa*gij2/(dens*dens*Cv*Cv);
  PetscReal tauE=1.0/sqrt(taut+taua+taud);

//  PetscPrintf(PETSC_COMM_WORLD,"taut %e \n",taut );
//  PetscPrintf(PETSC_COMM_WORLD,"taua %e \n",taua );
//  PetscPrintf(PETSC_COMM_WORLD,"tauM %e \n",tauM );
//  PetscPrintf(PETSC_COMM_WORLD,"tauE %e \n",tauE );


    tauC_bar = tauC;
    tauM_bar = tauM;
    tauE_bar = tauE;


   }





  //Transform to pressure primitive variables

  PetscReal TauTemp[4][4]={{0.0}};
  TauTemp[0][0] = tauC_bar;
  TauTemp[1][1] = tauM_bar;
  TauTemp[2][2] = tauM_bar;
  TauTemp[3][3] = tauE_bar;

  PetscInt aa,bb;
  for (aa=0;aa<dof;aa++)
 	   for (bb=0;bb<dof;bb++){
 		   tau[aa][bb] = A0inv[aa][0]*TauTemp[0][bb]
 		                +A0inv[aa][1]*TauTemp[1][bb]
 		                +A0inv[aa][2]*TauTemp[2][bb]
 		                +A0inv[aa][3]*TauTemp[3][bb];

 	  }
 //



  return 0;
}


//#undef __FUNCT__
//#define __FUNCT__ "Compute_inverse"
//PetscErrorCode Compute_inverse( PetscScalar (*Y)[4],PetscReal (*Yinv)[4])
//{
//
//
//
//    PetscReal DET;
//    PetscInt  i,j;
//    PetscInt  dof = 4;
//
//    DET =  Y[0][0]*(Y[1][1]*(Y[2][2]*Y[3][3]-Y[2][3]*Y[3][2])
//                   +Y[1][2]*(Y[2][3]*Y[3][1]-Y[2][1]*Y[3][3])
//                   +Y[1][3]*(Y[2][1]*Y[3][2]-Y[2][2]*Y[3][1]))
//          -Y[0][1]*(Y[1][0]*(Y[2][2]*Y[3][3]-Y[2][3]*Y[3][2])
//                   +Y[1][2]*(Y[2][3]*Y[3][0]-Y[2][0]*Y[3][3])
//                   +Y[1][3]*(Y[2][0]*Y[3][2]-Y[2][2]*Y[3][0]))
//          +Y[0][2]*(Y[1][0]*(Y[2][1]*Y[3][3]-Y[2][3]*Y[3][1])
//                   +Y[1][1]*(Y[2][3]*Y[3][0]-Y[2][0]*Y[3][3])
//                   +Y[1][3]*(Y[2][0]*Y[3][1]-Y[2][1]*Y[3][0]))
//          -Y[0][3]*(Y[1][0]*(Y[2][1]*Y[3][2]-Y[2][2]*Y[3][1])
//                   +Y[1][1]*(Y[2][2]*Y[3][0]-Y[2][0]*Y[3][2])
//                   +Y[1][2]*(Y[2][0]*Y[3][1]-Y[2][1]*Y[3][0]));
//
//
//   PetscReal COFACTOR[4][4]={{0.0}};
//
//
//
//    COFACTOR[0][0] = Y[1][1]*(Y[2][2]*Y[3][3]-Y[2][3]*Y[3][2])+Y[1][2]*(Y[2][3]*Y[3][1]-Y[2][1]*Y[3][3])+Y[1][3]*(Y[2][1]*Y[3][2]-Y[2][2]*Y[3][1]);
//    COFACTOR[0][1] = Y[1][0]*(Y[2][3]*Y[3][2]-Y[2][2]*Y[3][3])+Y[1][2]*(Y[2][0]*Y[3][3]-Y[2][3]*Y[3][0])+Y[1][3]*(Y[2][2]*Y[3][0]-Y[2][0]*Y[3][2]);
//    COFACTOR[0][2] = Y[1][0]*(Y[2][1]*Y[3][3]-Y[2][3]*Y[3][1])+Y[1][1]*(Y[2][3]*Y[3][0]-Y[2][0]*Y[3][3])+Y[1][3]*(Y[2][0]*Y[3][1]-Y[2][1]*Y[3][0]);
//    COFACTOR[0][3] = Y[1][0]*(Y[2][2]*Y[3][1]-Y[2][1]*Y[3][2])+Y[1][1]*(Y[2][0]*Y[3][2]-Y[2][2]*Y[3][0])+Y[1][2]*(Y[2][1]*Y[3][0]-Y[2][0]*Y[3][1]);
//    COFACTOR[1][0] = Y[0][1]*(Y[2][3]*Y[3][2]-Y[2][2]*Y[3][3])+Y[0][2]*(Y[2][1]*Y[3][3]-Y[2][3]*Y[3][1])+Y[0][3]*(Y[2][2]*Y[3][1]-Y[2][1]*Y[3][2]);
//    COFACTOR[1][1] = Y[0][0]*(Y[2][2]*Y[3][3]-Y[2][3]*Y[3][2])+Y[0][2]*(Y[2][3]*Y[3][0]-Y[2][0]*Y[3][3])+Y[0][3]*(Y[2][0]*Y[3][2]-Y[2][2]*Y[3][0]);
//    COFACTOR[1][2] = Y[0][0]*(Y[2][3]*Y[3][1]-Y[2][1]*Y[3][3])+Y[0][1]*(Y[2][0]*Y[3][3]-Y[2][3]*Y[3][0])+Y[0][3]*(Y[2][1]*Y[3][0]-Y[2][0]*Y[3][1]);
//    COFACTOR[1][3] = Y[0][0]*(Y[2][1]*Y[3][2]-Y[2][2]*Y[3][1])+Y[0][1]*(Y[2][2]*Y[3][0]-Y[2][0]*Y[3][2])+Y[0][2]*(Y[2][0]*Y[3][1]-Y[2][1]*Y[3][0]);
//    COFACTOR[2][0] = Y[0][1]*(Y[1][2]*Y[3][3]-Y[1][3]*Y[3][2])+Y[0][2]*(Y[1][3]*Y[3][1]-Y[1][1]*Y[3][3])+Y[0][3]*(Y[1][1]*Y[3][2]-Y[1][2]*Y[3][1]);
//    COFACTOR[2][1] = Y[0][0]*(Y[1][3]*Y[3][2]-Y[1][2]*Y[3][3])+Y[0][2]*(Y[1][0]*Y[3][3]-Y[1][3]*Y[3][0])+Y[0][3]*(Y[1][2]*Y[3][0]-Y[1][0]*Y[3][2]);
//    COFACTOR[2][2] = Y[0][0]*(Y[1][1]*Y[3][3]-Y[1][3]*Y[3][1])+Y[0][1]*(Y[1][3]*Y[3][0]-Y[1][0]*Y[3][3])+Y[0][3]*(Y[1][0]*Y[3][1]-Y[1][1]*Y[3][0]);
//    COFACTOR[2][3] = Y[0][0]*(Y[1][2]*Y[3][1]-Y[1][1]*Y[3][2])+Y[0][1]*(Y[1][0]*Y[3][2]-Y[1][2]*Y[3][0])+Y[0][2]*(Y[1][1]*Y[3][0]-Y[1][0]*Y[3][1]);
//    COFACTOR[3][0] = Y[0][1]*(Y[1][3]*Y[2][2]-Y[1][2]*Y[2][3])+Y[0][2]*(Y[1][1]*Y[2][3]-Y[1][3]*Y[2][1])+Y[0][3]*(Y[1][2]*Y[2][1]-Y[1][1]*Y[2][2]);
//    COFACTOR[3][1] = Y[0][0]*(Y[1][2]*Y[2][3]-Y[1][3]*Y[2][2])+Y[0][2]*(Y[1][3]*Y[2][0]-Y[1][0]*Y[2][3])+Y[0][3]*(Y[1][0]*Y[2][2]-Y[1][2]*Y[2][0]);
//    COFACTOR[3][2] = Y[0][0]*(Y[1][3]*Y[2][1]-Y[1][1]*Y[2][3])+Y[0][1]*(Y[1][0]*Y[2][3]-Y[1][3]*Y[2][0])+Y[0][3]*(Y[1][1]*Y[2][0]-Y[1][0]*Y[2][1]);
//    COFACTOR[3][3] = Y[0][0]*(Y[1][1]*Y[2][2]-Y[1][2]*Y[2][1])+Y[0][1]*(Y[1][2]*Y[2][0]-Y[1][0]*Y[2][2])+Y[0][2]*(Y[1][0]*Y[2][1]-Y[1][1]*Y[2][0]);
//
//
//
//
//            for (i=0;i<dof;i++){
//                for (j=0;j<dof;j++){
//                    Yinv[i][j]= COFACTOR[j][i]/DET;
//                }
//            }
//
//
//
//
//  return 0;
//}


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
  //PetscLogEventBegin(NS_Residual,0,0,0,0);
	PetscErrorCode ierr;
   AppCtx *user = (AppCtx *)ctx;


   PetscInt i,j,l;
   PetscInt  dof    = pnt->dof;
   PetscInt  dim    = pnt->dim;
   PetscReal lamda  = user->lamda;
   PetscReal mu  = user->mu;
   PetscReal kappa  = user->kappa;
   PetscReal Cv     = user->Cv;
   PetscReal RR     = user->R;
   PetscReal P0 = 100000;
   PetscReal B = 3000.0*P0;
   PetscReal N = 7.14;

   	 PetscScalar u[dof], u_t[dof];
     IGAPointFormValue(pnt,U,&u[0]);
     IGAPointFormValue(pnt,V,&u_t[0]);
     PetscScalar P   = u[0];
     PetscScalar ux  = u[1];
     PetscScalar uy  = u[2];
     PetscScalar temp= u[3];

    // PetscReal   mu = 0.0906*pow(temp,1.5)/(temp+0.0001406);
 //  PetscReal   s        = Cp*log(temp/temp0)-R*log(P/P0)+S0;
     PetscReal   alpha_p  = 1/temp;
     PetscReal   beta_t   = 1/P;
     PetscReal dens0    = 1000;
     PetscReal dens     = dens0*pow((1/B)*(P-P0)+1, 1/N);
     PetscReal f        = P0+B*(pow(dens/dens0,N)-1);
     PetscReal fprime   = (1.0/dens)*B*N*(pow(dens/dens0,N));
     PetscReal A011     = 1.0/fprime;
     PetscReal fcr      = P0+B*(pow(0.99,N)-1);
     if(P<=fcr+0.001){
       fprime = 0;
       dens = dens0*pow((1/B)*(fcr-P0)+1, 1/N);
       A011 = 0;
     }
     f        = P0+B*(pow(dens/dens0,N)-1);
     fprime   = (1.0/dens)*B*N*(pow(dens/dens0,N));
     A011     = 1.0/fprime;

     PetscReal   chi      = lamda+2*mu;

//   PetscScalar k    = (ux*ux+uy*uy+uz*uz)/2.0;
//   PetscScalar h    = Cp*temp;
//   PetscReal   e1   = h+k;
//   PetscReal   e1p  = dens*beta_t*e1 - alpha_p*temp;
//   PetscReal   e4p  =-dens*alpha_p*e1+dens*Cp;
//   PetscReal   e2p  = e1p+1.0; // Plus or minus??????
//   PetscReal   e3p  = dens*e1;
     PetscReal   kfac = 1/dens;
     PetscReal   ufac = ux*ux+uy*uy;
     PetscReal   etot = Cv*temp + 0.5*(ux*ux+uy*uy); // etot = cv*T + 1/2|u|^2


  PetscScalar grad_u[dof][dim];
  PetscScalar hess_u[dof][dim][dim];
  IGAPointFormGrad (pnt,U,&grad_u[0][0]);
  IGAPointFormHess(pnt,U,&hess_u[0][0][0]);
  PetscReal InvGradMap[dim][dim];
  IGAPointFormGradMap(pnt,NULL,&InvGradMap[0][0]);

  PetscReal        umi[2] = {0.0};
  PetscScalar   tau[4][4] = {{0.0}};
  PetscReal      A1_adv[4][4] = {{0.0}};
  PetscReal      A1_pt[4][4] = {{0.0}};
  PetscReal      A2_adv[4][4] = {{0.0}};
  PetscReal      A2_pt[4][4] = {{0.0}};
  PetscReal      A1_p[4][4] = {{0.0}};
  PetscReal      A2_p[4][4] = {{0.0}};
  PetscReal      A0[4][4] = {{0.0}};
  PetscReal   A0inv[4][4] = {{0.0}};
  PetscReal K[2][2][3][3] = {{{{0.0}}}};
  PetscReal       G[2][2] = {{0}};


    PetscReal t11,t12,t21,t22;
    t11 = mu*(grad_u[1][0]+grad_u[1][0])+lamda*(grad_u[1][0]+grad_u[2][1]);
    t22 = mu*(grad_u[2][1]+grad_u[2][1])+lamda*(grad_u[1][0]+grad_u[2][1]);
    t12 = mu*(grad_u[1][1]+grad_u[2][0]);
    t21 = mu*(grad_u[2][0]+grad_u[1][1]);



    A0[0][0] = A011;
    A0[1][0] = A011*ux;
    A0[1][1] = dens;
    A0[2][0] = A011*uy;
    A0[2][2] = dens;
    A0[3][0] = A011*Cv*temp;
    A0[3][3] = dens*Cv;


    A0inv[0][0] =  fprime;
    A0inv[1][0] = -ux/dens;
    A0inv[1][1] =  1.0/dens;
    A0inv[2][0] = -uy/dens;
    A0inv[2][2] =  1.0/dens;
    A0inv[3][0] =  -1.0/alpha_p/dens;
    A0inv[3][3] =  1.0/(Cv*dens);				//ok

    A1_adv[0][0] = A011*ux;
    A1_adv[0][1] = dens;
    A1_adv[1][0] = A011*ux*ux ;
    A1_adv[1][1] = 2.0*dens*ux;
    A1_adv[2][0] = A011*ux*uy;
    A1_adv[2][1] = dens*uy;
    A1_adv[2][2] = dens*ux;
    A1_adv[3][0] = ux*Cv*temp*A011;
    A1_adv[3][1] = dens*Cv*temp;
    A1_adv[3][3] = dens*ux*Cv*(-alpha_p*temp + 1.0);			//ok
    A1_p[1][0] = 1.0;
    A1_pt[3][1] = f-t11;
    A1_pt[3][2] = -t12;


  // Compute A2

  A2_adv[0][0] = A011*uy;
  A2_adv[0][2] = dens;
  A2_adv[1][0] = A011*ux*uy;
  A2_adv[1][1] = dens*uy;
  A2_adv[1][2] = dens*ux;
  A2_adv[2][0] = A011*uy*uy;
  A2_adv[2][2] = 2.0*dens*uy;
  A2_adv[3][0] = uy*Cv*temp*A011;
  A2_adv[3][2] = dens*Cv*temp;
  A2_adv[3][3] = dens*uy*Cv*(-alpha_p*temp + 1.0);	//ok
  A2_p[2][0] = 1.0;
  A2_pt[3][1]=-t21;
  A2_pt[3][2]=f-t22;



//  for (i=0;i<dof;i++)
//	for (j=0;j<dof;j++){
//		A1_c[i][j] -= umi[0]*A0[i][j];
//		A2_c[i][j] -= umi[1]*A0[i][j];
//	}



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
F1[1] = P;
F2[2] = P;



for (i=0;i<dim;i++)
for (j=0;j<dim;j++)
for (l=0;l<dim;l++){
G[i][j] += InvGradMap[i][l]*InvGradMap[j][l];
}

//ierr = Tau1(G,dt,u,tau,user,A0inv,umi);CHKERRQ(ierr);

PetscReal A1_cons[4][4] = {{0.0}};
PetscReal A2_cons[4][4] = {{0.0}};
PetscReal K_cons[2][2][4][4] = {{{{0.0}}}};

ComputeAMatrixConservation(u,A0inv,A1_adv,A1_p,A1_pt,A2_adv,A2_p,A2_pt,K,A1_cons,A2_cons,K_cons,user);

ierr = DirectTau(G,dt,u,tau,user,A0inv,A1_cons,A2_cons,K_cons,umi);CHKERRQ(ierr);

PetscInt m;
  //for (m=0;m<dof;m++){
   // for (j=0;j<dof;j++){
 // PetscPrintf(PETSC_COMM_WORLD,"tau %e \n",tau[m][j] );
   // }
 // }


  PetscReal  *N0 = pnt->shape[0];
  PetscReal (*N1)[dim] = (PetscReal (*)[dim]) pnt->shape[1];




  PetscReal Res[4] = {0.0};

  Res[0]  = -user->F[0];
  Res[1]  = -user->F[1];
  Res[2]  = -user->F[2];
  Res[3]  = -user->F[3];

  for (i=0;i<4;i++){
  for (j=0;j<4;j++){
  	    Res[i] += A0[i][j]*u_t[j] + (A1_adv[i][j]+A1_p[i][j]+A1_pt[i][j])*grad_u[j][0] + (A2_adv[i][j]+A2_p[i][j]+A2_pt[i][j])*grad_u[j][1];
  }
  }


  //Stabilization terms DC

   PetscReal hu[3] = {0.0};
   PetscReal A0gradY[4][2]={{0.0}};
   PetscReal tau_m,tau_c,tau_t;
   PetscReal Res_weighted_norm, GradU_weighted_norm, res_gradU_ratio, k_cap;

   Res_weighted_norm = fprime*sqrt(Res[0]*Res[0])+sqrt(ufac)*sqrt(Res[1]*Res[1]+Res[2]*Res[2])+sqrt(Res[3]+Res[3]);

   res_gradU_ratio = Res_weighted_norm;

   for (j=0;j<dim;j++)
   for (i=0;i<dof;i++)
   for (l=0;l<dof;l++)
   {
   	A0gradY[i][j] += A0[i][l]*grad_u[l][j];		//ok
   }


   //hu=0.0;
   for (i=0;i<dim;i++){
   for (j=0;j<dim;j++){
   	hu[0] += A0gradY[0][i]*G[i][j]*A0gradY[0][j];	//ok
   }
   }

   if (hu[0] > 1e-12){
   tau_c = sqrt(Res[0]*Res[0]/hu[0]);
   }else{
   tau_c = 0.0;
 }



   //hu=0.0;
   for (i=0;i<dim;i++){
   for (j=0;j<dim;j++){
   	hu[1] += A0gradY[1][i]*G[i][j]*A0gradY[1][j];   	//ok
   	hu[1] += A0gradY[2][i]*G[i][j]*A0gradY[2][j];		//ok
   }
   }

   if (hu[1] > 1e-12){
   tau_m = sqrt((Res[1]*Res[1] + Res[2]*Res[2]) /hu[1]);
   }else{
   tau_m = 0.0;
 }



   //hu=0.0;
   for (i=0;i<dim;i++)
   for (j=0;j<dim;j++){
   	hu[2] += A0gradY[3][i]*G[i][j]*A0gradY[3][j];	//ok
   }

   if (hu[2] > 1e-12){
   tau_t = sqrt(Res[3]*Res[3]/hu[2]);
   }else{
   tau_t = 0.0;
 }
//tau_c = res_gradU_ratio/(fprime*sqrt(hu[0])+sqrt(ufac)*sqrt(hu[1])+sqrt(hu[2])+1e-15);
//tau_m = tau_c;
//tau_t = tau_c;
   //PetscPrintf(PETSC_COMM_WORLD,"tau_c %e \n",tau_c );
   //PetscPrintf(PETSC_COMM_WORLD,"tau_m %e \n",tau_m );
   //PetscPrintf(PETSC_COMM_WORLD,"tau_t %e \n",tau_t );

   PetscReal Ginv[2][2]={{0.0}};
   PetscReal   det,detinv;

   det  = G[0][0]*G[1][1]-G[0][1]*G[1][0];

   detinv = 1.0/det;      //ok

   Ginv[0][0] =  G[1][1]*detinv;
   Ginv[1][0] = -G[1][0]*detinv;
   Ginv[0][1] = -G[0][1]*detinv;
   Ginv[1][1] =  G[0][0]*detinv; //ok

PetscReal uGu = 0.0;
for(int i = 0; i<dim; i++){
  for(int j =0; i<dim; i++){
    uGu += u[i+1]*Ginv[i][j]*u[j+1];
  }
}

k_cap = sqrt(fprime*(Ginv[0][0]+Ginv[1][1]+Ginv[2][2]));

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
    	   // R[a][i] += -N1[a][0]*A1_p[i][j]*u[j] - N1[a][1]*A2_p[i][j]*u[j];
		//R[a][i] += (A0[i][j]*u_t[j] + A1_c[i][j]*grad_u[j][0] + A2_c[i][j]*grad_u[j][1])*Na;
    }
    }

    for (i=0;i<4;i++){


    	    R[a][i] += -N1[a][0]*F1[i] - N1[a][1]*F2[i];

    }


//    for (m=0;m<dof;m++){
//    PetscPrintf(PETSC_COMM_WORLD,"Res1 %e \n",R[a][m] );
//    }

    for (l=0;l<dim;l++){
    for (m=0;m<dim;m++){
    	for (i=0;i<dof-1;i++){
        for (j=0;j<dof-1;j++){
            		R[a][i+1] += N1[a][l]*K[l][m][i][j]*grad_u[j+1][m];		//ok
        }
    	}
    }
    }


//    for (m=0;m<dof;m++){
//    PetscPrintf(PETSC_COMM_WORLD,"Res2 %e \n",R[a][m] );
//    }



    //Stabilization terms SUPG

// PetscReal RT[4] = {0.0};
    for (i=0;i<dof;i++)
    for (j=0;j<dof;j++)
    for (l=0;l<dof;l++)
    {
        R[a][i] += ((A1_adv[i][j]+A1_p[i][j]+A1_pt[i][j])*N1[a][0] + (A2_adv[i][j]+A2_p[i][j]+A2_pt[i][j])*N1[a][1])*tau[j][l]*Res[l];
	//RT[i] += ((A1_c[i][j]+A1_p[i][j])*N1[a][0] + (A2_c[i][j]+A2_p[i][j])*N1[a][1])*tau[j][l]*Res[l];
    }


   // for (m=0;m<dof;m++){
   // PetscPrintf(PETSC_COMM_WORLD,"RT %e \n",RT[m] );
   // }

//    PetscInt c,d;
//    for (i=0;i<dof;i++)
//               for (j=0;j<dof;j++)                 //Diffusion. Do we use it????
//             	  for (l=1;l<dof;l++)
//             	          for (m=1;m<dof;m++)
//                         	  for (c=0;c<dim;c++)
//                         	          for (d=0;d<dim;d++)  {
//             	        	  R[a][i] -= A1[i][j]*N1[a][0]*tau[j][l]*K[c][d][l-1][m-1]*hess_u[m][d][c];
//             	        	  R[a][i] -= A2[i][j]*N1[a][1]*tau[j][l]*K[c][d][l-1][m-1]*hess_u[m][d][c];
//  }






//    for (m=0;m<dof;m++){
//    PetscPrintf(PETSC_COMM_WORLD,"Res4 %e \n",R[a][m] );
//    }




//    PetscPrintf(PETSC_COMM_WORLD,"tau_c %e \n",tau_c );
//    PetscPrintf(PETSC_COMM_WORLD,"tau_t %e \n",tau_t );
//    PetscPrintf(PETSC_COMM_WORLD,"tau_t %e \n",tau_t );

//PetscReal RVDC[4]={0.0};
    PetscReal DC = 1.0;
    //if(DC*tau_c > k_cap || DC*tau_c != DC*tau_c){tau_c = k_cap;}
    //tau_m = tau_c;
    //tau_t = tau_c;

    for (i=0;i<dim;i++)
    {

    	    		  R[a][0] += DC*N1[a][i]*tau_c*A0gradY[0][i];
    	    		  R[a][1] += DC*N1[a][i]*tau_m*A0gradY[1][i];
    	    		  R[a][2] += DC*N1[a][i]*tau_m*A0gradY[2][i];
    	    		  R[a][3] += DC*N1[a][i]*tau_t*A0gradY[3][i];
    	    		 // RVDC[0] += N1[a][i]*tau_c*A0gradY[0][i];
    	    		  //RVDC[1] += N1[a][i]*tau_m*A0gradY[1][i];
    	    		 // RVDC[2] += N1[a][i]*tau_m*A0gradY[2][i];
    	    		 // RVDC[3] += N1[a][i]*tau_t*A0gradY[3][i];

    }



//for (m=0;m<dof;m++)    PetscPrintf(PETSC_COMM_WORLD,"RVDC %e \n",RVDC[m] );

//        for (m=0;m<dof;m++)    PetscPrintf(PETSC_COMM_WORLD,"uu %e \n",uu[m] );
//        for (m=0;m<dof;m++)    PetscPrintf(PETSC_COMM_WORLD,"vv %e \n",vv[m] );


//    for (m=0;m<5;m++){
//    for (j=0;j<5;j++){
//      PetscPrintf(PETSC_COMM_WORLD,"vdc %e \n",vdc[j][m] );
//    }
//    }
//    for (m=0;m<5;m++){
//      PetscPrintf(PETSC_COMM_WORLD,"N1 %e \n",N1[a][m] );
//    }
//
//    for (m=0;m<dim;m++){
//    for (j=0;j<dim;j++){
//    PetscPrintf(PETSC_COMM_WORLD,"G %e \n",G[m][j] );
//    }
//    }
//    for (m=0;m<dof;m++){
//    for (j=0;j<dim;j++){
//    PetscPrintf(PETSC_COMM_WORLD,"A0gradY %e \n",A0gradY[m][j] );
//    }
//    }

//
//      for (m=0;m<dof;m++){
//    //  for (j=0;j<dim;j++){
//      PetscPrintf(PETSC_COMM_WORLD,"Res5 %e \n",R[a][m] );
//    //  }
//      }


    }

  PetscInt n;
      if (pnt->parent->ID[0]==10){
      //	for (m=0;m<dof;m++){
      		//for (n=0;n<dof;n++){
      			//PetscPrintf(PETSC_COMM_WORLD,"tau %e\n",tau[m][n]);
//
//
      		//}
      	//}
      	//PetscPrintf(PETSC_COMM_WORLD,"tau_c %e \n",tau_c );
      //	PetscPrintf(PETSC_COMM_WORLD,"tau_m %e \n",tau_m );
      //	PetscPrintf(PETSC_COMM_WORLD,"tau_t %e \n",tau_t );
//
//
//
      	for (a=0; a<nen; a++){
      		for (i=0; i<dof; i++){
      			//PetscPrintf(PETSC_COMM_WORLD,"Residual %e \n",R[a][i] );

     		}
      	}
//
//
      		for (i=0; i<dof; i++){
//
     			//PetscPrintf(PETSC_COMM_WORLD,"Res %e \n",Res[i]);
      		}
//
      }


  //PetscLogEventEnd(NS_Residual,0,0,0,0);
  return 0;
}


#undef  __FUNCT__
#define __FUNCT__ "ResidualBound"
PetscErrorCode ResidualBound(IGAPoint pnt,PetscReal dt,
                        PetscReal shift,const PetscScalar *V,
                        PetscReal t,const PetscScalar *U,
                        PetscScalar *Re,void *ctx)
{
  //PetscLogEventBegin(NS_Residual,0,0,0,0);
	PetscErrorCode ierr;
   AppCtx *user = (AppCtx *)ctx;


   PetscInt i,j,l;
   PetscInt  dof    = pnt->dof;
   PetscInt  dim    = pnt->dim;
   PetscReal lamda  = user->lamda;

   PetscReal kappa  = user->kappa;
   PetscReal Cv     = user->Cv;
   PetscReal RR     = user->R;
   PetscReal P0 = 100000;
   PetscReal B = 3000.0*P0;
   PetscReal N = 7.14;

   	 PetscScalar u[dof], u_t[dof];
     IGAPointFormValue(pnt,U,&u[0]);
     IGAPointFormValue(pnt,V,&u_t[0]);
     PetscScalar P   = u[0];
     PetscScalar ux  = u[1];
     PetscScalar uy  = u[2];
     PetscScalar temp= u[3];


 //  PetscReal   s        = Cp*log(temp/temp0)-R*log(P/P0)+S0;
     PetscReal   alpha_p  = 1/temp;
     PetscReal   beta_t   = 1/P;
     PetscReal dens0    = 1000;
     PetscReal dens     = dens0*pow((1/B)*(P-P0)+1, 1/N);
     PetscReal f        = P0+B*(pow(dens/dens0,N)-1);
     PetscReal fprime   = (1.0/dens)*B*N*(pow(dens/dens0,N));
     PetscReal   A011     = 1.0/fprime;
     PetscReal fcr        = P0+B*(pow(0.99,N)-1);
     if(P<=fcr+0.001){
       fprime = 0;
       dens = dens0*pow((1/B)*(fcr-P0)+1, 1/N);
       A011 = 0;
     }
     f        = P0+B*(pow(dens/dens0,N)-1);
     fprime   = (1.0/dens)*B*N*(pow(dens/dens0,N));
     A011     = 1.0/fprime;
//   PetscScalar k    = (ux*ux+uy*uy+uz*uz)/2.0;
//   PetscScalar h    = Cp*temp;
//   PetscReal   e1   = h+k;
//   PetscReal   e1p  = dens*beta_t*e1 - alpha_p*temp;
//   PetscReal   e4p  =-dens*alpha_p*e1+dens*Cp;
//   PetscReal   e2p  = e1p+1.0; // Plus or minus??????
//   PetscReal   e3p  = dens*e1;
     PetscReal   kfac = 1/dens;
     PetscReal   ufac = ux*ux+uy*uy;
     PetscReal   etot = Cv*temp + 0.5*(ux*ux+uy*uy); // etot = cv*T + 1/2|u|^2


  PetscScalar grad_u[dof][dim];
  PetscScalar hess_u[dof][dim][dim];
  IGAPointFormGrad (pnt,U,&grad_u[0][0]);
  IGAPointFormHess(pnt,U,&hess_u[0][0][0]);
  PetscReal InvGradMap[dim][dim];
  IGAPointFormGradMap(pnt,NULL,&InvGradMap[0][0]);

  PetscReal        umi[2] = {0.0};

  PetscReal      A1_p[4][4] = {{0.0}};
  PetscReal      A2_p[4][4] = {{0.0}};
  PetscReal      A0[4][4] = {{0.0}};

  A0[0][0] = A011;
  A0[1][0] = A011*ux;
  A0[1][1] = dens;
  A0[2][0] = A011*uy;
  A0[2][2] = dens;
  A0[3][0] = A011*etot;
  A0[3][1] = dens*ux;
  A0[3][2] = dens*uy;
  A0[3][3] = dens*Cv;

	//ok

  A1_p[1][0] = 1.0;
  A1_p[3][0] = ux;
  A1_p[3][1] = f;


  A2_p[2][0] = 1.0;
  A2_p[3][0] = uy;
  A2_p[3][2] = f;


PetscReal F1[4]={0.0};
PetscReal F2[4]={0.0};
F1[1] = f;
F2[2] = f;


  //PetscPrintf(PETSC_COMM_WORLD,"element %d\n",pnt->parent->ID[0]);

  for (i=0;i<dof;i++)
	for (j=0;j<dof;j++){
		A1_p[i][j] -= umi[0]*A0[i][j];
		A2_p[i][j] -= umi[1]*A0[i][j];
	}

  PetscReal  *N0 = pnt->shape[0];
  PetscReal *norm;
  PetscInt m;

  norm = pnt->normal;

  //PetscPrintf(PETSC_COMM_WORLD,"norm1 %e\n",norm[0]);
// PetscPrintf(PETSC_COMM_WORLD,"norm2 %e\n",norm[1]);

  PetscScalar (*R)[dof] = (PetscScalar (*)[dof])Re;
  PetscInt a,nen=pnt->nen;
  for (a=0; a<nen; a++) {
    PetscReal Na    = N0[a];



    for (i=0;i<4;i++){
    //for (j=0;j<4;j++){
    	   // R[a][i] += Na*A1_p[i][j]*u[j]*norm[0] + Na*A2_p[i][j]*u[j]*norm[1];
		R[a][i] += Na*F1[i]*norm[0] + Na*F2[i]*norm[1];

    //}
    }


  }



 // if (pnt->parent->ID[0]==10){
	//for (a=0; a<nen; a++)
		//for (i=0; i<dof; i++)
			//PetscPrintf(PETSC_COMM_WORLD,"Resbound %e \n",R[a][i] );
 // }


  //PetscLogEventEnd(NS_Residual,0,0,0,0);
  return 0;
}




#undef  __FUNCT__
#define __FUNCT__ "Tangent"
PetscErrorCode Tangent(IGAPoint pnt,PetscReal dt,
                       PetscReal shift,const PetscScalar *V,
                       PetscReal t,const PetscScalar *U,
                       PetscScalar *Ke,void *ctx)
{
  //PetscLogEventBegin(NS_Tangent,0,0,0,0);

  AppCtx *user = (AppCtx *)ctx;


  PetscInt 	dof   = pnt->dof;
  PetscReal Cv    = user->Cv;
  PetscReal R     = user->R;
  PetscReal dim   = pnt->dim;
  PetscReal P0 = 100000;
  PetscReal B = 3000.0*P0;
  PetscReal N = 7.14;
//  PetscPrintf(PETSC_COMM_WORLD,"Tangent0 \n");

  PetscScalar u[dof];
  IGAPointFormValue(pnt,U,&u[0]);
  PetscReal *N0 = pnt->shape[0];



//  PetscPrintf(PETSC_COMM_WORLD,"P %e ux %e uy %e uz %e T %e\n",P,ux,uy,uz,temp);




  PetscReal A0[4][4]={{0.0}};


  PetscScalar P   = u[0];
  PetscScalar ux  = u[1];
  PetscScalar uy  = u[2];
  PetscScalar temp= u[3];


  PetscReal   etot = Cv*temp + 0.5*(ux*ux+uy*uy); // etot = cv*T + 1/2|u|^2
  PetscReal   alpha_p  = 1/temp;
  PetscReal   beta_t   = 1/P;
  PetscReal dens0    = 1000;
  PetscReal dens     = dens0*pow((1/B)*(P-P0)+1, 1/N);
  PetscReal f        = P0+B*(pow(dens/dens0,N)-1);
  PetscReal fprime   = (1.0/dens)*B*N*(pow(dens/dens0,N));
  PetscReal   A011     = 1.0/fprime;
  PetscReal fcr        = P0+B*(pow(0.99,N)-1);
  if(P<=fcr+0.001){
    fprime = 0.0;
    dens = dens0*pow((1/B)*(fcr-P0)+1, 1/N);
    A011 = 0.0;
  }
  f        = P0+B*(pow(dens/dens0,N)-1);
  fprime   = (1.0/dens)*B*N*(pow(dens/dens0,N));
  A011     = 1.0/fprime;

  PetscReal   kfac=1/dens;
  PetscReal   ufac=ux*ux+uy*uy;


	  PetscReal aa = dens*beta_t;
	  PetscReal bb =-dens*alpha_p;
	  PetscReal cc = dens*beta_t*etot;
  PetscReal ee = dens*(-alpha_p*etot + Cv);
  PetscReal factor = bb*cc-aa*ee;

  A0[0][0] = A011;
  A0[1][0] = A011*ux;
  A0[1][1] = dens;
  A0[2][0] = A011*uy;
  A0[2][2] = dens;
  A0[3][0] = A011*Cv*temp;
  A0[3][3] = dens*Cv;



			//ok

  PetscInt i,j;

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

//  PetscPrintf(PETSC_COMM_WORLD,"Tangent2 \n");

  //PetscLogEventEnd(NS_Tangent,0,0,0,0);
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

  PetscFunctionReturn(0);
}


//#undef  __FUNCT__
//#define __FUNCT__ "Tangent"
//PetscErrorCode Tangent(IGAPoint pnt,PetscReal dt,
//                       PetscReal shift,const PetscScalar *V,
//                       PetscReal t,const PetscScalar *U,
//                       PetscScalar *Ke,void *ctx)
//{
//  //PetscLogEventBegin(NS_Tangent,0,0,0,0);
//
//  AppCtx *user = (AppCtx *)ctx;
//
//
//  PetscInt 	dof   = pnt->dof;
//  PetscInt 	dim   = pnt->dim;
//  PetscReal lamda = user->lamda;
//  PetscReal mu    = user->mu;
//  PetscReal kappa = user->kappa;
//  PetscReal Cv    = user->Cv;
//  PetscReal Cp    = user->Cp;
//  PetscReal R     = user->R;
//
////  PetscPrintf(PETSC_COMM_WORLD,"Tangent0 \n");
//
//  PetscScalar u[dof];
//  IGAPointFormValue(pnt,U,&u[0]);
//  PetscScalar ux=u[1];
//  PetscScalar uy=u[2];
//  PetscScalar uz=u[3];
//  PetscScalar temp=u[4];
//
//
//  PetscScalar grad_u[dof][dim];
//  IGAPointFormGrad (pnt,U,&grad_u[0][0]);
//  PetscScalar ux_y=grad_u[1][1];
//  PetscScalar uy_z=grad_u[2][2];
//  PetscScalar uz_x=grad_u[3][0];
//
//
//  PetscReal   alpha_p  = 1/temp;
//  PetscReal   beta_t   = 1/u[0];
////  PetscReal   s        = Cp*log(temp/temp0)-R*log(u[0]/P0)+S0;
//  PetscScalar dens     = u[0]/(R*temp); ///???????
////  PetscScalar v        = 1.0/dens;
//
//  PetscScalar k   = (ux*ux+uy*uy+uz*uz)/2.0;
//  PetscScalar h   = Cp*temp;
//  PetscReal   e1  = h+k;
//  PetscReal   e1p = dens*beta_t*e1 - alpha_p*temp;
//  PetscReal   e4p =-dens*alpha_p*e1+dens*Cp;
//  PetscReal   e2p =e1p+1.0; // Plus or minus??????
//  PetscReal   e3p =dens*e1;
//  PetscReal   chi =lamda+2*mu;
//
//
//  PetscReal InvGradMap[dim][dim];
//  IGAPointFormGradMap(pnt,0,&InvGradMap[0][0]);
//
//   PetscReal A0_inv[dof][dof];
//         A0_inv[0][0]=(e4p+dens*alpha_p*(ux*ux+uy*uy+uz*uz))/(dens*beta_t*Cv);
//         A0_inv[1][0]=-ux;
//         A0_inv[2][0]=-uy;
//         A0_inv[3][0]=-uz;
//         A0_inv[4][0]=(dens*beta_t*(ux*ux+uy*uy+uz*uz)-e1p)/(dens*beta_t*Cv);
//
//         A0_inv[0][1]=-alpha_p*ux/(beta_t*Cv);
//         A0_inv[1][1]=1.0;
//         A0_inv[2][1]=0.0;
//         A0_inv[3][1]=0.0;
//         A0_inv[4][1]=-ux/Cv;
//
//         A0_inv[0][2]=-alpha_p*uy/(beta_t*Cv);
//         A0_inv[1][2]=0.0;
//         A0_inv[2][2]=1.0;
//         A0_inv[3][2]=0.0;
//         A0_inv[4][2]=-uy/Cv;
//
//         A0_inv[0][3]=-alpha_p*uz/(beta_t*Cv);
//         A0_inv[1][3]=0.0;
//         A0_inv[2][3]=0.0;
//         A0_inv[3][3]=1.0;
//         A0_inv[4][3]=-uz/Cv;
//
//         A0_inv[0][4]=alpha_p/(beta_t*Cv);
//         A0_inv[1][4]=0.0;
//         A0_inv[2][4]=0.0;
//         A0_inv[3][4]=0.0;
//         A0_inv[4][4]=1.0/Cv;
//
//         PetscInt i,j,l,m;
//
//  PetscReal *N0 = pnt->shape[0];
//  PetscReal (*N1)[dim] = (PetscReal (*)[dim]) pnt->shape[1];
//
//
//  PetscInt a,b,nen=pnt->nen;
//
//  PetscScalar (*KK)[dof][nen][dof] = (PetscScalar (*)[dof][nen][dof])Ke;
//  for (a=0; a<nen; a++) {
//
//
//      for (i=0;i<dof;i++)
//        for (j=0;j<dof;j++)
//        	KK[a][i][a][j] += A0_inv[i][j]/(dens*N0[a]);
//
////      for (i=0;i<dof;i++)
////              for (j=0;j<dof;j++)
////            	  for (l=0;l<dof;l++)
////            	          for (m=0;m<dof;m++){
////            	        	  KK[a][i][b][j] += A1[j][i]*N1[a][0]*tau[j][l]*A0[l][m]*N0[b];
////            	        	  KK[a][i][b][j] += A2[j][i]*N1[a][1]*tau[j][l]*A0[l][m]*N0[b];
////            	        	  KK[a][i][b][j] += A3[j][i]*N1[a][2]*tau[j][l]*A0[l][m]*N0[b];
////            	          }
//
//  }
//
////  PetscPrintf(PETSC_COMM_WORLD,"Tangent2 \n");
//
//  //PetscLogEventEnd(NS_Tangent,0,0,0,0);
//  return 0;
//}



typedef struct {
  PetscScalar p,ux,uy,temp;
} Field;



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
//      PetscReal x = (PetscReal)i / ( (PetscReal)(info.mx-1) )*user->Lx;
//      PetscReal y = (PetscReal)j / ( (PetscReal)(info.my-1) )*user->Ly;
        PetscReal x = iga->geometryX[((j-info.ys)*iga->geom_gwidth[0]+i-info.xs)*dim];
        PetscReal y = iga->geometryX[((j-info.ys)*iga->geom_gwidth[0]+i-info.xs)*dim + 1];
//        	PetscReal z = iga->geometryX[((j-info.zs)*iga->geom_gwidth[0]+i-info.zs)*dim + 1];

      //  PetscPrintf(PETSC_COMM_WORLD," %e\n",x);

//if (x <= 0.05){
//        u[j][i].ux   =  (-(y-0.5)*(y-0.5)+0.25)*4.0*0.0317;
//}else{
//	    u[j][i].ux   =  0;
//}
//        u[j][i].uy   =  0.0;
//        u[j][i].p    =  (1.5-0.5*x)*user->p0;
//        u[j][i].temp =  user->temp0;
//

       // if (x <= 2.1){
        //	u[j][i].ux   =  1.0;
         //  u[j][i].uy   =  0.0;
          //  u[j][i].p    =  user->p0;
          // u[j][i].temp =  user->temp0;
       // }else{
    	   // u[j][i].ux   =  0.375;
           // u[j][i].uy   =  0.0;
           // u[j][i].p    =  0.80357;
           // u[j][i].temp =  0.10453e-2;
       // }


//        u[j][i].ux   =  cos(10.0*3.14159/180.0);
//        u[j][i].uy   = -sin(10.0*3.14159/180.0);
//        u[j][i].p    = user->p0;
//        u[j][i].temp =  6.2306e-4;
//        if (y<=0.00001 && x>=0.00001)
//        	u[j][i].uy   =  0.0;
////



//        PetscReal r = sqrt(x*x + y*y);
//        PetscReal PI = 4.0*atan(1.0);
//
          PetscReal h = user->Lx/iga->elem_sizes[0];
//
//        u[j][i].ux = 0.0;
//        u[j][i].uy = 0.0;
//
//        if ((r>1.5556-h/2.0) && (r<1.5556+h/2.0)){
////        	u[j][i].temp = 1465.0;
////        	u[j][i].p = 14750000.00;
//        	u[j][i].temp =  -1.0/(user->Cv*h)*r + 1.0/user->Cv + 0.00002075;
//        	u[j][i].p = user->R*u[j][i].temp + 0.00595;
//        }else{
//        	u[j][i].temp = 0.00002075;
//        	u[j][i].p = 0.00595;
//        }


//
       // PetscInt index[1];
       // PetscReal value[1];
       // index[0] = 0.0;
       // VecGetValues(FluidPortion,1,index,value);

PetscReal xtemp = x+1.8;
PetscReal ytemp = y+1.8;
//PetscReal xtemp = x-3.6;
//PetscReal ytemp = y-3.7;
u[j][i].ux   =  -5.55136719*cos(ytemp/xtemp);
u[j][i].uy   =  -5.55136719*sin(ytemp/xtemp);
PetscReal f  = 100000+3000*100000*(pow(0.99,7.15)-1);
u[j][i].temp =  291;
u[j][i].ux   =  0.0;
u[j][i].uy   =  0.0;
u[j][i].p    =  f;
f  = 100000+3000*100000*(pow(1.001115,7.15)-1);
  if (xtemp*xtemp+ytemp*ytemp>=7.5*7.5)
	{
    u[j][i].ux   =  -5.55136719*cos(ytemp/xtemp);
    u[j][i].uy   =  -5.55136719*sin(ytemp/xtemp);
    u[j][i].p    =  f;
    //if(x==3.6){u[j][i].ux   =  -555.136719;}
    //if(y==3.6){u[j][i].uy   =  -555.136719;}
  }
  if (xtemp*xtemp+ytemp*ytemp>=7.5*7.5-0.01)
	{
    u[j][i].ux   =  -5.55136719*cos(ytemp/xtemp);
    u[j][i].uy   =  -5.55136719*sin(ytemp/xtemp);
  }
  if(u[j][i].ux!=u[j][i].ux) {exit(0);}
  if(u[j][i].uy!=u[j][i].uy) {exit(0);}



//  u[j][i].ux   =  0.0;
//  u[j][i].uy   =  0.0;
//  u[j][i].p    =  100000.0;
//  u[j][i].temp =  290.0;
////
//
//if (((y)*(y) + x*x) <= 0.001*0.001)
//{														// Kazem's example
//u[j][i].temp   = 1465.0;
//u[j][i].p   =    6750000.00;
//}




//  if (y<=0.00001 && x>=0.00001)
////    	  if (x>0.2){
//	u[j][i].uy   =  0.0;


//        u[j][i].ux   =  1195.9;
//        u[j][i].uy   =  0.0;
//        u[j][i].p    =  100000;
//        u[j][i].temp =  291;
//
//
//           if   ((x<=0.014) && (x>=0.003) && (y<=0.364*x+0.0139) && (y>=-0.364*x+0.01609)){
//        	      u[j][i].ux   =  0.0;
//        	      u[j][i].uy   =  0.0;
//           }


//           u[j][i].ux   =  1195.9;
//           u[j][i].uy   =  0.0;
//           u[j][i].p    =  100000;                                                              //Fixed cylinder
//           u[j][i].temp =  291;
//
//              if   ((x - 0.01)*(x - 0.01) + (y - 0.015)*(y - 0.015) <= 0.0025*0.0025){
//                      u[j][i].ux   =  0.0;
//                      u[j][i].uy   =  0.0;
//              }








//        u[j][i].ux   =  0.0;
//              u[j][i].uy   =  0.0;
//              u[j][i].p    =  100000.0;
//              u[j][i].temp =  290.0;
//  //
//
//        if (((y-0.2)*(y-0.2) + x*x) <= 0.061*0.061)
//      	{
//      	u[j][i].temp   = 1465.0;
//          u[j][i].p   =    6746268.65;
//        }



//    	u[j][i].ux   =  1.0;
//        u[j][i].uy   = 0.0;
//        u[j][i].p    =  0.079365079;
//        u[j][i].temp =  2.769e-4;
//
//  if (y<=0.00001)
//	  if (x>0.2){
//  	u[j][i].ux   =  0.0;
//    u[j][i].temp =  7.754e-4;
//  }




//        PetscPrintf(PETSC_COMM_WORLD," %d, %d\n",i,j);


    }
  }
  ierr = DMDAVecRestoreArray(da,U,&u);CHKERRQ(ierr);
  ierr = DMDestroy(&da);;CHKERRQ(ierr);
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
static PetscErrorCode TSPredictStage_GeneralizedAlpha(AppCtx *user)
{

//  PetscPrintf(PETSC_COMM_WORLD,"					TSPredictStage_JB\n");

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
static PetscErrorCode TSUpdateAlphaLevels_GeneralizedAlpha(AppCtx *user)
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


	   Vec cont;
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


	   //VecSet(dens,0.0);

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

	 //	  PetscPrintf(PETSC_COMM_WORLD,"nodes: %d \n",Nnodes);
	 //	  PetscPrintf(PETSC_COMM_WORLD,"Nloc: %d \n",Nloc);
	 //	  PetscPrintf(PETSC_COMM_WORLD,"dof: %d \n",dof);
	 //	  PetscPrintf(PETSC_COMM_WORLD,"Nloc/dof: %d \n",Nloc/dof);
	 //
	 //	  PetscInt NlocD,NlocV, NlocT;
	 //	  ierr  = VecGetLocalSize(dens,&NlocD);CHKERRQ(ierr);
	 //	  ierr  = VecGetLocalSize(vel,&NlocV);CHKERRQ(ierr);
	 //	  ierr  = VecGetLocalSize(temp,&NlocT);CHKERRQ(ierr);
	 //	  PetscPrintf(PETSC_COMM_WORLD,"NlocD: %d \n",NlocD);
	 //	  PetscPrintf(PETSC_COMM_WORLD,"NlocV: %d \n",NlocV);
	 //	  PetscPrintf(PETSC_COMM_WORLD,"NlocT: %d \n",NlocT);
	 //	  PetscInt NglobD,NglobV, NglobT;
	 //	  ierr  = VecGetSize(dens,&NglobD);CHKERRQ(ierr);
	 //	  ierr  = VecGetSize(vel,&NglobV);CHKERRQ(ierr);
	 //	  ierr  = VecGetSize(temp,&NglobT);CHKERRQ(ierr);
	 //	  PetscPrintf(PETSC_COMM_WORLD,"NglobD: %d \n",NglobD);
	 //	  PetscPrintf(PETSC_COMM_WORLD,"NglobV: %d \n",NglobV);
	 //	  PetscPrintf(PETSC_COMM_WORLD,"NglobT: %d \n",NglobT);


	   ierr = VecNormBegin(cont, NORM_2, &normCont);CHKERRQ(ierr);
	   ierr = VecNormEnd(cont, NORM_2, &normCont);CHKERRQ(ierr);

	   ierr = VecNormBegin(mom, NORM_2, &normMom);CHKERRQ(ierr);
	   ierr = VecNormEnd(mom, NORM_2, &normMom);CHKERRQ(ierr);

	   ierr = VecNormBegin(ener, NORM_2, &normEner);CHKERRQ(ierr);
	   ierr = VecNormEnd(ener, NORM_2, &normEner);CHKERRQ(ierr);

//      PetscPrintf(PETSC_COMM_WORLD,"it: %d, Res norm Cont: %e  \n",it,normCont);
//      PetscPrintf(PETSC_COMM_WORLD,"it: %d, Res norm Mom: %e  \n",it,normMom);
//      PetscPrintf(PETSC_COMM_WORLD,"it: %d, Res norm Ener: %e  \n",it,normEner);


      ierr = VecDestroy(&cont);CHKERRQ(ierr);
      ierr = VecDestroy(&mom);CHKERRQ(ierr);
      ierr = VecDestroy(&ener);CHKERRQ(ierr);



  PetscFunctionReturn(0);
}


//#undef __FUNCT__
//#define __FUNCT__ "ComputeResNorm"
//static PetscErrorCode ComputeResNorm(IGA iga,Vec Res, PetscInt it, PetscReal t, PetscInt step, PetscInt maxit)
//{
//
//
//  PetscErrorCode ierr;
//
//  PetscFunctionBegin;
//  PetscInt        dof = iga->dof;
//  PetscInt        dim = iga->dim;
//  PetscReal           fnorm;
//  PetscInt i,j,Nloc;
//  VecGetLocalSize(Res,&Nloc);
//  PetscInt Nnodes = iga->node_sizes[0]*iga->node_sizes[1]*iga->node_sizes[2];
//
//
//
//      //Global Residual norm (including all dof)
//
//      ierr = VecNormBegin(Res,NORM_2,&fnorm);CHKERRQ(ierr);        // fnorm <- ||F||
//      ierr = VecNormEnd(Res,NORM_2,&fnorm);CHKERRQ(ierr);
////      PetscPrintf(PETSC_COMM_WORLD,"it: %d, Res norm: %e  \n",it,fnorm);
//
//
//
//
//
//
//
//
//      //Residual norm by dof
//
//
//       PetscReal           normCont;
//       PetscReal           normMom;
//       PetscReal           normEner;
//
//
//       Vec cont;
//       VecCreate(PETSC_COMM_WORLD,&cont);
//       VecSetSizes(cont,Nloc/dof,Nnodes);
//       VecSetFromOptions(cont);
//
//       Vec mom;
//       VecCreate(PETSC_COMM_WORLD,&mom);
//       VecSetSizes(mom,Nloc/dof*dim,Nnodes*dim);
//       VecSetFromOptions(mom);
//
//       Vec ener;
//       VecCreate(PETSC_COMM_WORLD,&ener);
//       VecSetSizes(ener,Nloc/dof,Nnodes);
//       VecSetFromOptions(ener);
//
//
//       //VecSet(dens,0.0);
//
//       PetscScalar *arrayRes;
//       PetscScalar *aCont;
//       PetscScalar *aMom;
//       PetscScalar *aEner;
//
//       ierr  = VecGetArray(Res,&arrayRes);CHKERRQ(ierr);
//
//       ierr  = VecGetArray(cont,&aCont);CHKERRQ(ierr);
//       ierr  = VecGetArray(mom,&aMom);CHKERRQ(ierr);
//       ierr  = VecGetArray(ener,&aEner);CHKERRQ(ierr);
//
//
//          for(i=0;i<Nloc/dof;i++) {
//            aCont[i]    = arrayRes[dof*i];
//            aEner[i]    = arrayRes[dof*i + dof - 1];
//               for(j=0;j<dim;j++){
//                 aMom[dim*i + j]    = arrayRes[dof*i + 1 + j];
//               }
//          }
//
//
//      ierr  = VecRestoreArray(ener,&aEner);CHKERRQ(ierr);
//      ierr  = VecRestoreArray(mom,&aMom);CHKERRQ(ierr);
//      ierr  = VecRestoreArray(cont,&aCont);CHKERRQ(ierr);
//      ierr  = VecRestoreArray(Res,&arrayRes);CHKERRQ(ierr);
//
//     //      PetscPrintf(PETSC_COMM_WORLD,"nodes: %d \n",Nnodes);
//     //      PetscPrintf(PETSC_COMM_WORLD,"Nloc: %d \n",Nloc);
//     //      PetscPrintf(PETSC_COMM_WORLD,"dof: %d \n",dof);
//     //      PetscPrintf(PETSC_COMM_WORLD,"Nloc/dof: %d \n",Nloc/dof);
//     //
//     //      PetscInt NlocD,NlocV, NlocT;
//     //      ierr  = VecGetLocalSize(dens,&NlocD);CHKERRQ(ierr);
//     //      ierr  = VecGetLocalSize(vel,&NlocV);CHKERRQ(ierr);
//     //      ierr  = VecGetLocalSize(temp,&NlocT);CHKERRQ(ierr);
//     //      PetscPrintf(PETSC_COMM_WORLD,"NlocD: %d \n",NlocD);
//     //      PetscPrintf(PETSC_COMM_WORLD,"NlocV: %d \n",NlocV);
//     //      PetscPrintf(PETSC_COMM_WORLD,"NlocT: %d \n",NlocT);
//     //      PetscInt NglobD,NglobV, NglobT;
//     //      ierr  = VecGetSize(dens,&NglobD);CHKERRQ(ierr);
//     //      ierr  = VecGetSize(vel,&NglobV);CHKERRQ(ierr);
//     //      ierr  = VecGetSize(temp,&NglobT);CHKERRQ(ierr);
//     //      PetscPrintf(PETSC_COMM_WORLD,"NglobD: %d \n",NglobD);
//     //      PetscPrintf(PETSC_COMM_WORLD,"NglobV: %d \n",NglobV);
//     //      PetscPrintf(PETSC_COMM_WORLD,"NglobT: %d \n",NglobT);
//
//
//       ierr = VecNormBegin(cont, NORM_2, &normCont);CHKERRQ(ierr);
//       ierr = VecNormEnd(cont, NORM_2, &normCont);CHKERRQ(ierr);
//
//       ierr = VecNormBegin(mom, NORM_2, &normMom);CHKERRQ(ierr);
//       ierr = VecNormEnd(mom, NORM_2, &normMom);CHKERRQ(ierr);
//
//       ierr = VecNormBegin(ener, NORM_2, &normEner);CHKERRQ(ierr);
//       ierr = VecNormEnd(ener, NORM_2, &normEner);CHKERRQ(ierr);
//
////      PetscPrintf(PETSC_COMM_WORLD,"it: %d, Res norm Cont: %e  \n",it,normCont);
////      PetscPrintf(PETSC_COMM_WORLD,"it: %d, Res norm Mom: %e  \n",it,normMom);
////      PetscPrintf(PETSC_COMM_WORLD,"it: %d, Res norm Ener: %e  \n",it,normEner);
//
//
//
////          char filenameEnergy[256];
////           sprintf(filenameEnergy,"Norm.txt");
////
////           FILE *fL = fopen(filenameEnergy, "a");
////           if (fL == NULL)
////             {
////             printf("Error opening file!\n");
////             exit(1);
////             }
////           if(step == 0 && it==0)   PetscFPrintf(PETSC_COMM_WORLD,fL,"Time, Total ResNorm, ResNormP, ResNormV, ResNormE \n");
////           if(it==maxit-1)        PetscFPrintf(PETSC_COMM_WORLD,fL,"%e, %e, %e, %e, %e \n", t, fnorm, normCont, normMom, normEner);
////               fclose(fL);
//
//
//      ierr = VecDestroy(&cont);CHKERRQ(ierr);
//      ierr = VecDestroy(&mom);CHKERRQ(ierr);
//      ierr = VecDestroy(&ener);CHKERRQ(ierr);
//
//
//
//  PetscFunctionReturn(0);
//}




#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char *argv[]) {


  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);

  AppCtx user;



  user.mu       =  0.0;//3.137e-3;//1.983e-5;
  user.lamda    =  0.0;
  user.kappa    = 0.0; //0.024;
  user.Cv       = 4137.9;
  user.R        = 286.6;
  user.Cp       = 4181.1;

 //user.temp0    = 0.61941e-3;//0.25;
 //user.p0       = 0.177; //71.0;
//  user.Lx    = 1.0;   user.Ly    = 1.0;

  //user.temp0    = 6.2306e-4;//0.25;
  user.p0       = 100000; //71.0;

  user.Lx    = 3.6;   user.Ly    = 3.6;



  user.F[0]=0.0;
  user.F[1]=0.0;
  user.F[2]=0.0;
  user.F[3]=0.0;


  user.currentTime  =0.0;
  user.timeStep		=1.0e-5;
  user.finalTime 	=0.006;
  user.stepNumber 	=0;

  user.max_its = 3;


  user.FreqResults = 1;
  user.TimeRestart = 0;
  user.StepRestart = 0;
  user.FreqRestarts= 1000;

//  PetscBool output = PETSC_FALSE;
//  PetscBool monitor = PETSC_FALSE;


  IGA iga;
  ierr = IGACreate(PETSC_COMM_WORLD,&iga);CHKERRQ(ierr);
  ierr = IGASetDim(iga,2);CHKERRQ(ierr);
  ierr = IGASetDof(iga,4);CHKERRQ(ierr);
  user.iga = iga;

//  ierr = IGARead(iga,"./Geo/Immersed_square_30x30_1.0x1.0.dat");CHKERRQ(ierr);
  ierr = IGARead(iga,"./Geo/Geometry.dat");CHKERRQ(ierr);

//  IGABoundary bnd;
//  PetscInt dir=1,side,field;
//  for (side=0;side<2;side++) {
//    ierr = IGAGetBoundary(iga,dir,side,&bnd);CHKERRQ(ierr);
//    for (field=0;field<5;field++) {
//      ierr = IGABoundarySetValue(bnd,field,100.0);CHKERRQ(ierr);
//    }
//  }

  ierr = IGASetFromOptions(iga);CHKERRQ(ierr);
  ierr = IGASetUp(iga);CHKERRQ(ierr);
  ierr = IGAWrite(iga,"igaF.dat");CHKERRQ(ierr);



  //Create solution vector V (velocities) and A (accelerations)
  PetscReal t=0;
  ierr = IGACreateVec(iga,&user.V0);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&user.A0);CHKERRQ(ierr);
  ierr = VecZeroEntries(user.V0);CHKERRQ(ierr);
  ierr = VecZeroEntries(user.A0);CHKERRQ(ierr);
  ierr = FormInitialCondition(iga,t,user.V0,&user);CHKERRQ(ierr);

  ierr = IGACreateVec(iga,&user.dA);CHKERRQ(ierr);
  ierr = VecZeroEntries(user.dA);CHKERRQ(ierr);




  // Dump Initial Solution
  char filename[256];
  sprintf(filename,"velS%d.dat",user.stepNumber);
  ierr = IGAWriteVec(user.iga,user.V0,filename);CHKERRQ(ierr);




//  MPI_Comm comm;
//  PetscViewer viewer;
//  ierr = PetscObjectGetComm((PetscObject)V,&comm);CHKERRQ(ierr);
//  char filenameRestart[256];
//  sprintf(filenameRestart,"Solution%d.dat",user.stepNumber);
//  ierr = PetscViewerBinaryOpen(comm,filenameRestart,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
//  ierr = VecView(V,viewer);CHKERRQ(ierr);
//  ierr = PetscViewerDestroy(&viewer);


  IGABoundary bound;



  ierr = IGASetUserIFunction(iga,Residual,&user);CHKERRQ(ierr);
  ierr =IGAGetBoundary(iga,0,0,&bound);CHKERRQ(ierr);
  ierr =IGABoundarySetUserIFunction(bound,ResidualBound,&user);CHKERRQ(ierr);
  ierr =IGAGetBoundary(iga,0,1,&bound);CHKERRQ(ierr);
  ierr =IGABoundarySetUserIFunction(bound,ResidualBound,&user);CHKERRQ(ierr);
  ierr =IGAGetBoundary(iga,1,0,&bound);CHKERRQ(ierr);
  ierr =IGABoundarySetUserIFunction(bound,ResidualBound,&user);CHKERRQ(ierr);
  ierr =IGAGetBoundary(iga,1,1,&bound);CHKERRQ(ierr);
  ierr =IGABoundarySetUserIFunction(bound,ResidualBound,&user);CHKERRQ(ierr);
  ierr = IGASetUserIJacobian(iga,Tangent,&user);CHKERRQ(ierr);


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

  //######################################
  //         	TIME STEP
  //######################################

  while(user.currentTime+user.timeStep <= user.finalTime+1e-9)
   {


	  	 PetscPrintf(PETSC_COMM_WORLD,"######################################################### \n");
	  	 PetscPrintf(PETSC_COMM_WORLD,"Step Number: %d  Time step: %e, Time: %e \n",user.stepNumber, user.timeStep, user.currentTime);
	  	 PetscPrintf(PETSC_COMM_WORLD,"######################################################### \n");



	      ierr = TSPredictStage_GeneralizedAlpha(&user);CHKERRQ(ierr);


	         PetscInt maxits = user.max_its;               /* maximum number of iterations */

	         PetscInt i,j,k,m,l,it;
	         PetscInt dof = iga->dof;
	         PetscInt dim = iga->dim;
	         Mat A0;
	         Vec Res,SP;
	         ierr = IGACreateMat(iga,&A0);CHKERRQ(ierr);
	         ierr = IGACreateVec(iga,&Res);CHKERRQ(ierr);
	         ierr = IGACreateVec(iga,&SP);CHKERRQ(ierr);  //Could we compute SP in the first TS and stored?



	         for (it=0; it<maxits; it++) {




	       	  PetscPrintf(PETSC_COMM_WORLD,":::::::::::::::::::::::::::::::::\n");
	       	  PetscPrintf(PETSC_COMM_WORLD,"  Iteration: %d  \n", it);
	       	  PetscPrintf(PETSC_COMM_WORLD,":::::::::::::::::::::::::::::::::\n");



	       	  ierr = TSUpdateAlphaLevels_GeneralizedAlpha(&user);CHKERRQ(ierr);


//	       		           ierr = VecView(user.Va,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

	             PetscReal t = user.currentTime;
	             PetscReal dt = user.timeStep;
	             PetscReal Alpha_m = user.Alpha_m;
	             PetscReal Alpha_f = user.Alpha_f;
	             PetscReal stage_time = t + Alpha_f*dt;
	       //      th->scale_F = 1; //
	       //      th->scale_J = Alpha_f*Beta*dt*dt;
	       //      th->shift_V = Gamma/(dt*Beta);
	       //      th->shift_A = Alpha_m/(Alpha_f*dt*dt*Beta);


	           PetscPrintf(PETSC_COMM_WORLD,"IGAComputeIFunction \n");
	           ierr = IGAComputeIFunction(iga,dt,Alpha_m,user.Aa,stage_time,user.Va,Res);CHKERRQ(ierr);
	           PetscPrintf(PETSC_COMM_WORLD,"IGAComputeIJacobianComp \n");
	           //Do we need to compute SP in every it and time step. Should we stored??
	           ierr = IGAComputeIJacobianComp(iga,dt,Alpha_m,user.Aa,stage_time,user.Va,A0,&user);CHKERRQ(ierr);









	           PetscInt numFluidNodes  = iga->geom_lwidth[0]*iga->geom_lwidth[1]*iga->geom_lwidth[2];
    	        PetscInt nodesX  = iga->geom_lwidth[0], nodesY  = iga->geom_lwidth[1], nodesZ  = iga->geom_lwidth[2];
    	        PetscInt gnodesX = iga->geom_gwidth[0], gnodesY = iga->geom_gwidth[1];



      	         for(m=0;m<nodesZ;m++) {
      	             for(l=0;l<nodesY;l++) {
      	                 for(k=0;k<nodesX;k++) {

      	           /*if(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim + 0] <=0.00001) {
      	            	PetscInt index = (m*gnodesX*gnodesY + l*gnodesX+ k)*dof+1;
      	            	PetscInt index_array[4]={0.0};
      	            	index_array[0]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof;
      	            	index_array[1]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+1;
      	            	index_array[2]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+2;
      	            	index_array[3]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+3;

      	            	MatSetValuesLocal(A0,1,&index,4,index_array,bc_ux,INSERT_VALUES);CHKERRQ(ierr);
      	            	MatSetValuesLocal(A0,4,index_array,1,&index,bc_ux,INSERT_VALUES);CHKERRQ(ierr);
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


//      	            if(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim + 0] <=0.00001) {
//      	            	PetscInt index = (m*gnodesX*gnodesY + l*gnodesX+ k)*dof+2;
//      	            	PetscInt index_array[4]={0.0};
//      	            	index_array[0]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof;
//      	            	index_array[1]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+1;
//      	            	index_array[2]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+2;
//      	            	index_array[3]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+3;
//
//      	            	MatSetValuesLocal(A0,1,&index,4,index_array,bc_uy,INSERT_VALUES);CHKERRQ(ierr);
//      	            	MatSetValuesLocal(A0,4,index_array,1,&index,bc_uy,INSERT_VALUES);CHKERRQ(ierr);
//      	            	VecSetValueLocal(Res,index,0.0,INSERT_VALUES);CHKERRQ(ierr);
//
//      	            }
//
//
//      	            if(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim + 0] <=0.00001) {
//      	            	PetscInt index = (m*gnodesX*gnodesY + l*gnodesX+ k)*dof+3;
//      	            	PetscInt index_array[4]={0.0};
//      	            	index_array[0]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof;
//      	            	index_array[1]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+1;
//      	            	index_array[2]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+2;
//      	            	index_array[3]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+3;
//
//      	            	MatSetValuesLocal(A0,1,&index,4,index_array,bc_temp,INSERT_VALUES);CHKERRQ(ierr);
//      	            	MatSetValuesLocal(A0,4,index_array,1,&index,bc_temp,INSERT_VALUES);CHKERRQ(ierr);
//      	            	VecSetValueLocal(Res,index,0.0,INSERT_VALUES);CHKERRQ(ierr);
//
//      	            }

      	            /*if(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim +1] <=0.00001) {
      	            	PetscInt index = (m*gnodesX*gnodesY + l*gnodesX+ k)*dof+1;
      	            	PetscInt index_array[4]={0.0};
      	            	index_array[0]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof;
      	            	index_array[1]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+1;
      	            	index_array[2]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+2;
      	            	index_array[3]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+3;

      	            	MatSetValuesLocal(A0,1,&index,4,index_array,bc_ux,INSERT_VALUES);CHKERRQ(ierr);
      	            	MatSetValuesLocal(A0,4,index_array,1,&index,bc_ux,INSERT_VALUES);CHKERRQ(ierr);
      	            	VecSetValueLocal(Res,index,0.0,INSERT_VALUES);CHKERRQ(ierr);

      	            }*/

//      	            if((iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim] >=0.003) && (iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim] <=0.014) && (iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim+1] <=0.364*iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim]+0.0139) && (iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim+1] >=-0.364*iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim]+0.01609)) {//(user.Lx -0.00001)) {
//      	             	            	PetscInt index = (m*gnodesX*gnodesY + l*gnodesX+ k)*dof+1;
//      	             	            	PetscInt index_array[4]={0.0};
//      	             	            	index_array[0]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof;
//      	             	            	index_array[1]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+1;
//      	             	            	index_array[2]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+2;
//      	             	            	index_array[3]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+3;
//
//      	             	            	MatSetValuesLocal(A0,1,&index,4,index_array,bc_ux,INSERT_VALUES);CHKERRQ(ierr);
//      	             	            	MatSetValuesLocal(A0,4,index_array,1,&index,bc_ux,INSERT_VALUES);CHKERRQ(ierr);
//      	             	            	VecSetValueLocal(Res,index,0.0,INSERT_VALUES);CHKERRQ(ierr);
//
//      	             	            }
//      	            if((iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim] >=0.003) && (iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim] <=0.014) && (iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim+1] <=0.364*iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim]+0.0139) && (iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim+1] >=-0.364*iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim]+0.01609)) {//(user.Lx -0.00001)) {
//      	             	            	PetscInt index = (m*gnodesX*gnodesY + l*gnodesX+ k)*dof+2;
//      	             	            	PetscInt index_array[4]={0.0};
//      	             	            	index_array[0]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof;
//      	             	            	index_array[1]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+1;
//      	             	            	index_array[2]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+2;
//      	             	            	index_array[3]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+3;
//
//      	             	            	MatSetValuesLocal(A0,1,&index,4,index_array,bc_uy,INSERT_VALUES);CHKERRQ(ierr);
//      	             	            	MatSetValuesLocal(A0,4,index_array,1,&index,bc_uy,INSERT_VALUES);CHKERRQ(ierr);
//      	             	            	VecSetValueLocal(Res,index,0.0,INSERT_VALUES);CHKERRQ(ierr);
//
//      	             	            }



   //                           if((iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim] - 0.01)*(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim] - 0.01) +  (iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim+1] - 0.015)*(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim+1] - 0.015) <= 0.0025*0.0025) {
   //
   //
   //                                      PetscInt index = (m*gnodesX*gnodesY + l*gnodesX+ k)*dof+1;
   //                                      PetscInt index_array[4]={0.0};
   //                                      index_array[0]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof;
   //                                      index_array[1]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+1;
   //                                      index_array[2]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+2;
   //                                      index_array[3]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+3;
   //
   //                                      MatSetValuesLocal(A0,1,&index,4,index_array,bc_ux,INSERT_VALUES);CHKERRQ(ierr);
   //                                      MatSetValuesLocal(A0,4,index_array,1,&index,bc_ux,INSERT_VALUES);CHKERRQ(ierr);
   //                                      VecSetValueLocal(Res,index,0.0,INSERT_VALUES);CHKERRQ(ierr);
   //
   //                           }
   //
   //
   //
   //
   //                           if((iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim] - 0.01)*(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim] - 0.01) + (iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim+1] - 0.015)*(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim+1] - 0.015) <= 0.0025*0.0025) {
   //
   //
   //                                      PetscInt index = (m*gnodesX*gnodesY + l*gnodesX+ k)*dof+2;
   //                                      PetscInt index_array[4]={0.0};
   //                                      index_array[0]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof;
   //                                      index_array[1]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+1;
   //                                      index_array[2]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+2;
   //                                      index_array[3]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+3;
   //
   //                                      MatSetValuesLocal(A0,1,&index,4,index_array,bc_uy,INSERT_VALUES);CHKERRQ(ierr);
   //                                      MatSetValuesLocal(A0,4,index_array,1,&index,bc_uy,INSERT_VALUES);CHKERRQ(ierr);
   //                                      VecSetValueLocal(Res,index,0.0,INSERT_VALUES);CHKERRQ(ierr);
   //
   //                           }



//       	            if((iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim + 1] <=0.00001) &&  (iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim] >=0.2)){
//       	            	PetscInt index = (m*gnodesX*gnodesY + l*gnodesX+ k)*dof+1;
//       	            	PetscInt index_array[4]={0.0};
//       	            	index_array[0]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof;
//       	            	index_array[1]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+1;
//       	            	index_array[2]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+2;
//       	            	index_array[3]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+3;
//
//       	            	MatSetValuesLocal(A0,1,&index,4,index_array,bc_ux,INSERT_VALUES);CHKERRQ(ierr);
//       	            	MatSetValuesLocal(A0,4,index_array,1,&index,bc_ux,INSERT_VALUES);CHKERRQ(ierr);
//       	            	VecSetValueLocal(Res,index,0.0,INSERT_VALUES);CHKERRQ(ierr);
//
//       	            }
//
//       	            if((iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim + 1] <=0.00001) &&  (iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim] >=0.2)){
//       	            	PetscInt index = (m*gnodesX*gnodesY + l*gnodesX+ k)*dof+3;
//       	            	PetscInt index_array[4]={0.0};
//       	            	index_array[0]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof;
//       	            	index_array[1]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+1;
//       	            	index_array[2]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+2;
//       	            	index_array[3]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+3;
//
//       	            	MatSetValuesLocal(A0,1,&index,4,index_array,bc_temp,INSERT_VALUES);CHKERRQ(ierr);
//       	            	MatSetValuesLocal(A0,4,index_array,1,&index,bc_temp,INSERT_VALUES);CHKERRQ(ierr);
//       	            	VecSetValueLocal(Res,index,0.0,INSERT_VALUES);Cif(P<HKERRQ(ierr);
//
//       	            }
//
       	            /*if(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim +1] >=(user.Ly -0.00001)) {
       	            	PetscInt index = (m*gnodesX*nodesY + l*gnodesX+ k)*dof+1;
       	            	PetscInt index_array[4]={0.0};
       	            	index_array[0]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof;
      	            	index_array[1]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+1;
       	            	index_array[2]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+2;
       	            	index_array[3]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+3;

       	            	MatSetValuesLocal(A0,1,&index,4,index_array,bc_ux,INSERT_VALUES);CHKERRQ(ierr);
      	            	MatSetValuesLocal(A0,4,index_array,1,&index,bc_ux,INSERT_VALUES);CHKERRQ(ierr);
     	            	VecSetValueLocal(Res,index,0.0,INSERT_VALUES);CHKERRQ(ierr);

                  }*/
//   ////
//       	            if(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim +1] >=(user.Ly -0.00001)) {
//       	            	PetscInt index = (m*gnodesX*gnodesY + l*gnodesX+ k)*dof+1;
//       	            	PetscInt index_array[4]={0.0};
//       	            	index_array[0]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof;
//       	            	index_array[1]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+1;
//       	            	index_array[2]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+2;
//       	            	index_array[3]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+3;
//
//       	            	MatSetValuesLocal(A0,1,&index,4,index_array,bc_ux,INSERT_VALUES);CHKERRQ(ierr);
//       	            	MatSetValuesLocal(A0,4,index_array,1,&index,bc_ux,INSERT_VALUES);CHKERRQ(ierr);
//       	            	VecSetValueLocal(Res,index,0.0,INSERT_VALUES);CHKERRQ(ierr);
//
//       	            }
      	            /*if(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim +1] >=(user.Ly -0.00001)) {
      	            	PetscInt index = (m*gnodesX*gnodesY + l*gnodesX+ k)*dof+2;
      	            	PetscInt index_array[4]={0.0};
      	            	index_array[0]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof;
      	            	index_array[1]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+1;
      	            	index_array[2]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+2;
      	            	index_array[3]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+3;

      	            	MatSetValuesLocal(A0,1,&index,4,index_array,bc_uy,INSERT_VALUES);CHKERRQ(ierr);
      	            	MatSetValuesLocal(A0,4,index_array,1,&index,bc_uy,INSERT_VALUES);CHKERRQ(ierr);
      	            	VecSetValueLocal(Res,index,0.0,INSERT_VALUES);CHKERRQ(ierr);

      	            }*/
//       	            if(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim +1] >=(user.Ly -0.00001)) {
//       	            	PetscInt index = (m*gnodesX*gnodesY + l*gnodesX+ k)*dof+3;
//       	            	PetscInt index_array[4]={0.0};
//       	            	index_array[0]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof;
//       	            	index_array[1]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+1;
//       	            	index_array[2]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+2;
//       	            	index_array[3]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+3;
//
//       	            	MatSetValuesLocal(A0,1,&index,4,index_array,bc_temp,INSERT_VALUES);CHKERRQ(ierr);
//       	            	MatSetValuesLocal(A0,4,index_array,1,&index,bc_temp,INSERT_VALUES);CHKERRQ(ierr);
//       	            	VecSetValueLocal(Res,index,0.0,INSERT_VALUES);CHKERRQ(ierr);
//
//       	            }
//
     	            /*if(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim] >=(user.Lx -0.00001)) {
      	             	            	PetscInt index = (m*gnodesX*gnodesY + l*gnodesX+ k)*dof+1;
      	             	            	PetscInt index_array[4]={0.0};
      	             	            	index_array[0]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof;
      	             	            	index_array[1]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+1;
      	             	            	index_array[2]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+2;
      	             	            	index_array[3]=(m*gnodesX*gnodesY + l*gnodesX+ k)*dof+3;

      	             	            	MatSetValuesLocal(A0,1,&index,4,index_array,bc_ux,INSERT_VALUES);CHKERRQ(ierr);
      	             	            	MatSetValuesLocal(A0,4,index_array,1,&index,bc_ux,INSERT_VALUES);CHKERRQ(ierr);
      	             	            	VecSetValueLocal(Res,index,0.0,INSERT_VALUES);CHKERRQ(ierr);

      	             	            }*/



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


//	           ierr = VecView(SP,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	       //    ierr = MatView(A0inv,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);







	                 	 	    //Boundary Conditions
//	                 	 	   	ierr  = VecGetArray(user.dA,&arraydA);CHKERRQ(ierr);



//
//	                 	         for(m=0;m<nodesZ;m++) {
//	                 	             for(l=0;l<nodesY;l++) {
//	                 	                 for(k=0;k<nodesX;k++) {
//
//	                 	              		if(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim + 0] <=0.00001)  arraydA[(m*nodesX*nodesY + l*nodesX+ k)*dof + 0] = 0.0;
//	                 	             		if(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim + 0] <=0.00001)  arraydA[(m*nodesX*nodesY + l*nodesX+ k)*dof + 1] = 0.0;
//	                 	             		if(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim + 0] <=0.00001)  arraydA[(m*nodesX*nodesY + l*nodesX+ k)*dof + 2] = 0.0;
//	                 	             		if(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim + 0] <=0.00001)  arraydA[(m*nodesX*nodesY + l*nodesX+ k)*dof + 3] = 0.0;
//
//	                 	             		//if(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim + 0] >= (user.Lx -0.00001))  arraydA[(m*nodesX*nodesY + l*nodesX+ k)*dof + 0] = 0.0;
//	                 	             		//if(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim + 0] >= (user.Lx -0.00001))  arraydA[(m*nodesX*nodesY + l*nodesX+ k)*dof + 1] = 0.0;
//	                 	             		//if(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim + 0] >= (user.Lx -0.00001))  arraydA[(m*nodesX*nodesY + l*nodesX+ k)*dof + 2] = 0.0;
//	                 	             		//if(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim + 0] >= (user.Lx -0.00001))  arraydA[(m*nodesX*nodesY + l*nodesX+ k)*dof + 3] = 0.0;
//
////	                 	             		if(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim + 1] <=0.00001)  arraydA[(m*nodesX*nodesY + l*nodesX+ k)*dof + 0] = 0.0;
////	                 	             		if(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim + 1] <=0.00001)  arraydA[(m*nodesX*nodesY + l*nodesX+ k)*dof + 1] = 0.0;
//	                 	             		if(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim + 1] <=0.00001)  arraydA[(m*nodesX*nodesY + l*nodesX+ k)*dof + 2] = 0.0;
////	                 	             		if(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim + 1] <=0.00001)  arraydA[(m*nodesX*nodesY + l*nodesX+ k)*dof + 3] = 0.0;
//
//	                 	             		if(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim + 1] >= (user.Ly -0.00001))  arraydA[(m*nodesX*nodesY + l*nodesX+ k)*dof + 0] = 0.0;
//	                 	             		if(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim + 1] >= (user.Ly -0.00001))  arraydA[(m*nodesX*nodesY + l*nodesX+ k)*dof + 1] = 0.0;
//                 	             	    	if(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim + 1] >= (user.Ly -0.00001))  arraydA[(m*nodesX*nodesY + l*nodesX+ k)*dof + 2] = 0.0;
//	                 	             		if(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim + 1] >= (user.Ly -0.00001))  arraydA[(m*nodesX*nodesY + l*nodesX+ k)*dof + 3] = 0.0;
////

//if(iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim + 1] <=0.00001)
//	if (iga->geometryX[(m*gnodesX*gnodesY + l*gnodesX+k)*dim + 0] >0.2){
//	arraydA[(m*nodesX*nodesY + l*nodesX+ k)*dof + 1] = 0.0;
//	arraydA[(m*nodesX*nodesY + l*nodesX+ k)*dof + 3] = 0.0;
//}
//
//	                 	                 }
//	                 	             }
//	                 	         }



//		             	           PetscPrintf(PETSC_COMM_WORLD,"dA \n");
//		             	           ierr = VecView(user.dA,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
//		             	           PetscPrintf(PETSC_COMM_WORLD,"V1 \n");
//		             	           ierr = VecView(user.V1,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

		                 	 	   PetscReal Ynorm;
		                 	 	  ierr = VecNormBegin(user.dA,NORM_2,&Ynorm);CHKERRQ(ierr);        // fnorm <- ||F||
		                 	 	  ierr = VecNormEnd(user.dA,NORM_2,&Ynorm);CHKERRQ(ierr);
		                 	 	  PetscPrintf(PETSC_COMM_WORLD,"it: %d, Y norm: %e  \n",it,Ynorm);



	           ierr = TSUpdateStage_GeneralizedAlpha(&user,user.dA);CHKERRQ(ierr);


	           ierr = KSPDestroy(&ksp);CHKERRQ(ierr);

//

	         }  //End iteration loop



	           ierr = VecCopy(user.V1,user.V0);CHKERRQ(ierr);
	           ierr = VecCopy(user.A1,user.A0);CHKERRQ(ierr);

//	         ierr = VecView(user.V1,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

	    user.currentTime+=user.timeStep;
	    user.stepNumber++;


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





//    MPI_Comm comm;
//    PetscViewer viewer;
//    ierr = PetscObjectGetComm((PetscObject)V,&comm);CHKERRQ(ierr);
//    char filenameRestart[256];
//    sprintf(filenameRestart,"Solution%d.dat",user.stepNumber);
//    ierr = PetscViewerBinaryOpen(comm,filenameRestart,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
//    ierr = VecView(V,viewer);CHKERRQ(ierr);
//    ierr = PetscViewerDestroy(&viewer);






    ierr = MatDestroy(&A0);CHKERRQ(ierr);
    ierr = VecDestroy(&Res);CHKERRQ(ierr);
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







   ierr = IGADestroy(&iga);CHKERRQ(ierr);



  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}
