/****************************** turbulence model header file *******************************/
/***************************************** chi86 *******************************************/
/**************************************** Oct 2019 *****************************************/
/*******************************************************************************************/
#include "TurbulenceModel.h"

/*
 * muT & aT based on the k-epsilon-v2-f-model
 */
void kEv2f(Mesh *mesh, scalar CmuKv2f_0,
	               scalar sigmaE_0,
	               scalar CE2_0,
	               scalar CL_0,
	               scalar Ceta_0,
	               scalar C1_0,
	               scalar C2_0) {
  
  /* Medic, G. and Durbin, P.A., */
  /* "Towards improved prediction of heat transfer on turbine blades", */
  /* ASME, J. Turbomach. 2012. */
  
  int i;

  scalar RxK,RxE,RxV2;
  scalar CN;

  scalar CmuKv2f,sigmaK,sigmaE,sigmaV2,CE2,CT,CL,Ceta,C1,C2;

  scalar N;

  scalar boundkL,boundkR,boundeL,boundV2L,boundV2R,boundfL;

  /* scalar t,dt; */

  field l,m,r,rhs;//,x;

  scalar it;

  scalar epsilonWall,fWall;

  scalar muFp,muFm,tauFp,tauFm;
  
  scalar PV2Exp,CrossV2;
  scalar PEExp,CrossE;
  scalar PKExp,CrossK;


  /* scalar residuum; */


  /* residuum=0.0; */

  
  l   = mesh->l;
  m   = mesh->m;
  r   = mesh->r;
  rhs = mesh->rhs;
  
  //x = mesh->x;  
  

  RxK  =0.6;
  RxE  =0.6;
  RxV2 =0.6;

  CN=0.; // fully implicit

  sigmaK  = 1.0;
  sigmaV2 = 1.0;
  CT      = 6.0;

  // defined as "model constant"
  /* CmuKv2f = 0.22; */
  /* sigmaE  = 1.3; */
  /* CE2     = 1.9; */
  /* CL      = 0.23; */
  /* Ceta    = 70.0; */
  /* C1      = 1.4; */
  /* C2      = 0.3; */

  CmuKv2f = CmuKv2f_0;
  sigmaE  = sigmaE_0;
  CE2     = CE2_0;
  CL      = CL_0;
  Ceta    = Ceta_0;
  C1      = C1_0;
  C2      = C2_0;


  

  N=6.0;
  
  // boundary conditions
  boundkL = 1.0;                                          // Neumann @ center
  boundkR = - mesh->cells[imax+1].dx/mesh->cells[imax].dx;  // Dirichlet @ wall

  boundeL = 1.0;                                          // Neumann @ center

  boundV2L = 1.0;                                         // Neumann @ center
  boundV2R = - mesh->cells[imax+1].dx/mesh->cells[imax].dx; //  Dirichlet @ wall

  boundfL = 1.0;                                          // Neumann @ center


  /* t=0.0; */
  /* dt=1; */




  // bounding
  epsilonWall=mesh->EWall*mesh->cells[imax].nu*mesh->cells[imax].k;
  epsilonWall=fmin(epsilonWall,100.0);
 
    
        
  fWall=-4.0*(6.0-N)*powf(mesh->cells[imax].nu,2.0)*mesh->cells[imax].V2/(mesh->cells[imax].epsilon*powf(mesh->cells[imax].y,4.0))/powf(ReTau,2.0);

  //printf("%f %f %f\n",epsilonWall,fWall,mesh->cells[imax].k);

  

  mesh->cells[0].k=mesh->cells[1].k;
  mesh->cells[imax+1].k=-mesh->cells[imax].k;

  mesh->cells[0].epsilon=mesh->cells[1].epsilon;
  mesh->cells[imax].epsilon=epsilonWall;
  mesh->cells[imax+1].epsilon=epsilonWall;

  mesh->cells[0].V2=mesh->cells[1].V2;
  mesh->cells[imax+1].V2=-mesh->cells[imax].V2;

  mesh->cells[0].f=mesh->cells[1].f;
  mesh->cells[imax].f=fWall;
  mesh->cells[imax+1].f=fWall;





  mesh->cells[0].strainR=(mesh->cells[1].w-mesh->cells[0].w)/(mesh->cells[1].xc-mesh->cells[0].xc);
  for(i=1;i<=imax;i+=1) {
    mesh->cells[i].strainR=(mesh->cells[i+1].w-mesh->cells[i-1].w)/(mesh->cells[i+1].xc-mesh->cells[i-1].xc);
  }
  mesh->cells[imax+1].strainR=(mesh->cells[imax+1].w-mesh->cells[imax].w)/(mesh->cells[imax+1].xc-mesh->cells[imax].xc);


  mesh->cells[0].muT=PREC_MIN;
  mesh->cells[0].aT=PREC_MIN;
  for(i=0;i<=imax+1;i+=1) {
    mesh->cells[i].CE1=1.4*(1.0+0.05*sqrt(mesh->cells[i].k/mesh->cells[i].V2));
    
    mesh->cells[i].T=fmax(mesh->cells[i].k/mesh->cells[i].epsilon,	\
    			 CT*sqrt(mesh->cells[i].nu/(mesh->cells[i].epsilon*ReTau)));
    
    /* mesh->cells[i].T=fmin( fmax(mesh->cells[i].k/mesh->cells[i].epsilon, \ */
    /* 				CT*sqrt(mesh->cells[i].nu/(mesh->cells[i].epsilon*ReTau))), \ */
    /* 			   0.6*mesh->cells[i].k/(sqrt(6.0)*mesh->cells[i].V2*CmuKv2f*abs(mesh->cells[i].strainR))); */

    
    mesh->cells[i].L=CL*fmax(powf(mesh->cells[i].k,3.0/2.0)/mesh->cells[i].epsilon, \
    			    Ceta*powf(mesh->cells[i].nu,3.0/4.0)/(powf(mesh->cells[i].epsilon,1.0/4.0)*powf(ReTau,3.0/4.0)));
    
    /* mesh->cells[i].L=CL*fmax( fmin(powf(mesh->cells[i].k,3.0/2.0)/mesh->cells[i].epsilon, \ */
    /* 				   powf(mesh->cells[i].k,3.0/2.0)/(sqrt(6.0)*mesh->cells[i].V2*CmuKv2f*abs(mesh->cells[i].strainR))), \ */
    /* 			      Ceta*powf(mesh->cells[i].nu,3.0/4.0)/(powf(mesh->cells[i].epsilon,1.0/4.0)*powf(ReTau,3.0/4.0))); */
    
    mesh->cells[i].muT=fmax(mesh->cells[i].rho*CmuKv2f*mesh->cells[i].V2*mesh->cells[i].T*ReTau, PREC_MIN);
    mesh->cells[i].aT=mesh->cells[i].muT/mesh->cells[i].PrT;

    mesh->cells[i].PExpTurb=mesh->cells[i].muT*powf(mesh->cells[i].strainR,2.0);
  }
  mesh->cells[imax+1].muT=mesh->cells[imax].muT;
  mesh->cells[imax+1].aT=mesh->cells[imax].aT;




  
  // f -loop
  for(i=1;i<=imax;i+=1) {
    l[i-1]= - mesh->cells[i].xFm * powf(mesh->cells[i].L,2.0)/(mesh->cells[i].xc-mesh->cells[i-1].xc);
    
    m[i-1]= 1.0 + powf(mesh->cells[i].L,2.0)*				\
      (mesh->cells[i].xFp /(mesh->cells[i+1].xc - mesh->cells[i  ].xc)  + \
       mesh->cells[i].xFm /(mesh->cells[i  ].xc - mesh->cells[i-1].xc));

    r[i-1]= - mesh->cells[i].xFp * powf(mesh->cells[i].L,2.0) /(mesh->cells[i+1].xc-mesh->cells[i].xc);

    rhs[i-1]=(C1-1.0)*(2.0/3.0 - mesh->cells[i].V2/mesh->cells[i].k)/mesh->cells[i].T + C2*mesh->cells[i].PExpTurb/(mesh->cells[i].k*mesh->cells[i].rho*ReTau) + (N-1.0)*mesh->cells[i].V2/(mesh->cells[i].k*mesh->cells[i].T);
  }
  
  ThomasAlgWallModel(imax-1,l,m,r,rhs,boundfL,fWall);

  
  for(i=1;i<=imax;i+=1) {
    it=rhs[i-1];
    mesh->cells[i].f=fmax(it,PREC_MIN);
  }

 
  //  v2 -loop
  for(i=1;i<=imax;i+=1) {
    muFp=( (mesh->cells[i+1].mu+mesh->cells[i  ].mu)*0.5 + (mesh->cells[i+1].muT+mesh->cells[i  ].muT)*0.5/sigmaV2 )/ ReTau;
    muFm=( (mesh->cells[i  ].mu+mesh->cells[i-1].mu)*0.5 + (mesh->cells[i  ].muT+mesh->cells[i-1].muT)*0.5/sigmaV2 )/ ReTau;

    tauFp=muFp*(mesh->cells[i+1].V2-mesh->cells[i  ].V2)/(mesh->cells[i+1].xc-mesh->cells[i  ].xc);
    tauFm=muFm*(mesh->cells[i  ].V2-mesh->cells[i-1].V2)/(mesh->cells[i  ].xc-mesh->cells[i-1].xc);

    PV2Exp = mesh->cells[i].rho*mesh->cells[i].f*mesh->cells[i].k;
    CrossV2= - N*mesh->cells[i].rho*mesh->cells[i].epsilon/mesh->cells[i].k;

    l[i-1]= -(1.0-CN) * mesh->cells[i].xFm * muFm/(mesh->cells[i].xc-mesh->cells[i-1].xc);
    
    m[i-1]= ddt*mesh->cells[i].rho +					\
      (1.0-CN) *							\
      (mesh->cells[i].xFp * muFp/(mesh->cells[i+1].xc - mesh->cells[i  ].xc)  + \
       mesh->cells[i].xFm * muFm/(mesh->cells[i  ].xc - mesh->cells[i-1].xc)) - \
      (1.0-CN)*CrossV2;

    r[i-1]= -(1.0-CN) * mesh->cells[i].xFp *muFp/(mesh->cells[i+1].xc-mesh->cells[i].xc);

    //x[i-1]=mesh->cells[i  ].V2;

    rhs[i-1]=CN*CrossV2*mesh->cells[i].V2 + PV2Exp +			\
      CN*(mesh->cells[i].xFp*tauFp - mesh->cells[i].xFm*tauFm) +	\
      mesh->cells[i].rho*mesh->cells[i].V2*ddt;
  }
  
  ThomasAlg(imax-1,l,m,r,rhs,boundV2L,boundV2R);
  //ThomasAlgRelax(imax-1,l,m,r,rhs,x,boundV2L,boundV2R,RxV2);

  for(i=1;i<=imax;i+=1) {
    it=RxV2*rhs[i-1]+(1-RxV2)*mesh->cells[i].rho*mesh->cells[i].V2;
    //it=rhs[i-1];
    mesh->cells[i].V2=fmax(it/mesh->cells[i].rho,PREC_MIN);
  }


  //  epsilon -loop
  for(i=1;i<=imax;i+=1) {
    muFp=( (mesh->cells[i+1].mu+mesh->cells[i  ].mu)*0.5 + (mesh->cells[i+1].muT+mesh->cells[i  ].muT)*0.5/sigmaE )/ ReTau;
    muFm=( (mesh->cells[i  ].mu+mesh->cells[i-1].mu)*0.5 + (mesh->cells[i  ].muT+mesh->cells[i-1].muT)*0.5/sigmaE )/ ReTau;

    tauFp=muFp*(mesh->cells[i+1].epsilon-mesh->cells[i  ].epsilon)/(mesh->cells[i+1].xc-mesh->cells[i  ].xc);
    tauFm=muFm*(mesh->cells[i  ].epsilon-mesh->cells[i-1].epsilon)/(mesh->cells[i  ].xc-mesh->cells[i-1].xc);

    PEExp  = mesh->cells[i].CE1*mesh->cells[i].PExpTurb*mesh->cells[i].rho/(mesh->cells[i].T*ReTau);
    CrossE = - CE2*mesh->cells[i].rho/mesh->cells[i].T;
    

    l[i-1]= -(1.0-CN) * mesh->cells[i].xFm * muFm/(mesh->cells[i].xc-mesh->cells[i-1].xc);
    
    m[i-1]= ddt*mesh->cells[i].rho +					\
      (1.0-CN) *							\
      (mesh->cells[i].xFp * muFp/(mesh->cells[i+1].xc - mesh->cells[i  ].xc)  + \
       mesh->cells[i].xFm * muFm/(mesh->cells[i  ].xc - mesh->cells[i-1].xc)) -	\
      (1.0-CN)*CrossE;

    r[i-1]= -(1.0-CN) * mesh->cells[i].xFp *muFp/(mesh->cells[i+1].xc-mesh->cells[i].xc);

    //x[i-1]=mesh->cells[i  ].epsilon;

    rhs[i-1]=CN*CrossE*mesh->cells[i].epsilon + PEExp +			\
      CN*(mesh->cells[i].xFp*tauFp - mesh->cells[i].xFm*tauFm) +	\
      mesh->cells[i].rho*mesh->cells[i].epsilon*ddt;

    
  } 
  ThomasAlgWallModel(imax-1,l,m,r,rhs,boundeL,epsilonWall);
  //ThomasAlgWallModelRelax(imax-1,l,m,r,rhs,x,boundeL,epsilonWall,RxE);

  for(i=1;i<=imax;i+=1) {
    it=RxE*rhs[i-1]+(1.0-RxE)*mesh->cells[i].rho*mesh->cells[i].epsilon;
    //it=rhs[i-1];
    mesh->cells[i].epsilon=fmax(it/mesh->cells[i].rho,PREC_MIN);

  }




  //  k -loop
  for(i=1;i<=imax;i+=1) {
    muFp=( (mesh->cells[i+1].mu+mesh->cells[i  ].mu)*0.5 + (mesh->cells[i+1].muT+mesh->cells[i  ].muT)*0.5/sigmaK )/ ReTau;
    muFm=( (mesh->cells[i  ].mu+mesh->cells[i-1].mu)*0.5 + (mesh->cells[i  ].muT+mesh->cells[i-1].muT)*0.5/sigmaK )/ ReTau;

    tauFp=muFp*(mesh->cells[i+1].k-mesh->cells[i  ].k)/(mesh->cells[i+1].xc-mesh->cells[i  ].xc);
    tauFm=muFm*(mesh->cells[i  ].k-mesh->cells[i-1].k)/(mesh->cells[i  ].xc-mesh->cells[i-1].xc);

    PKExp = mesh->cells[i].PExpTurb/ReTau;
    CrossK   =  - mesh->cells[i].rho*mesh->cells[i].epsilon/mesh->cells[i].k;

    l[i-1]= -(1.0-CN) * mesh->cells[i].xFm * muFm/(mesh->cells[i].xc-mesh->cells[i-1].xc);
    
    m[i-1]= ddt*mesh->cells[i].rho +				\
      (1.0-CN) *							\
      (mesh->cells[i].xFp * muFp/(mesh->cells[i+1].xc - mesh->cells[i  ].xc)  + \
       mesh->cells[i].xFm * muFm/(mesh->cells[i  ].xc - mesh->cells[i-1].xc)) - \
      (1.0-CN)*CrossK;

    r[i-1]= -(1.0-CN) * mesh->cells[i].xFp *muFp/(mesh->cells[i+1].xc-mesh->cells[i].xc);

    //x[i-1]=mesh->cells[i  ].k;
    
    rhs[i-1]=CN*CrossK*mesh->cells[i].k + PKExp +			\
      CN*(mesh->cells[i].xFp*tauFp - mesh->cells[i].xFm*tauFm) +	\
      mesh->cells[i].rho*mesh->cells[i].k*ddt;

    
  }
  ThomasAlg(imax-1,l,m,r,rhs,boundkL,boundkR);
  //ThomasAlgRelax(imax-1,l,m,r,rhs,x,boundkL,boundkR,RxK);

  for(i=1;i<=imax;i+=1) {
    it=RxK*rhs[i-1]+(1.0-RxK)*mesh->cells[i].rho*mesh->cells[i].k;
    //it=rhs[i-1];
    mesh->cells[i].k=fmax(it/mesh->cells[i].rho,PREC_MIN);
    mesh->cells[i].Zeta=mesh->cells[i].V2/mesh->cells[i].k;
  }

  
  //printf("\t k/epsilon: residuum= %10.5e\n",residuum);


  /* mesh->cells[0].muT=1E-12; */
  /* for(i=1;i<=imax;i+=1) { */
  /*   mesh->cells[i].T=fmax(mesh->cells[i].k/mesh->cells[i].epsilon,	\ */
  /* 			 CT*sqrt(mesh->cells[i].nu/(mesh->cells[i].epsilon*ReTau))); */
  /*   mesh->cells[i].muT=fmax(mesh->cells[i].rho*CmuKv2f*mesh->cells[i].V2*mesh->cells[i].T*ReTau, PREC_MIN); */
  /* } */
  /* mesh->cells[imax+1].muT=mesh->cells[imax].muT; */
  
}
