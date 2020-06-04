/*                                                                *
 *                    Various CFD related functions               *
 *                                                                *
 *                                                                *
 *                    chi86                                       *
 *                    12.02.2020                                  *
 *                                                                *
 */



#include "CFD.h"
#include "Modules.h"

/*
 *  Initialize flow files and scalars
 */
void InitializeFlow(Mesh *mesh) {
  int i;
  int SW,SWH;
  
  SW=1;
  SWH=1;
  
  for(i=0;i<=imax+1;i+=1) {
    if( Wlog(kappa,beta,fabs(mesh->cells[i].yp)) < Wlin(fabs(mesh->cells[i].yp)) && SW==1 ) {
      mesh->cells[i].w=Wlog(kappa,beta,fabs(mesh->cells[i].yp));
    }
    else {
      SW=0;
      mesh->cells[i].w=Wlin(fabs(mesh->cells[i].yp));
    }
  
    if( (Hlog(kappa,beta,PrW,PrT0,fabs(mesh->cells[i].yp))<Hlin(PrW,fabs(mesh->cells[i].yp))) && SWH==1) {
      mesh->cells[i].h=Hlog(kappa,beta,PrW,PrT0,fabs(mesh->cells[i].yp));
    }
    else {
      SWH=0;
      mesh->cells[i].h=Hlin(PrW,fabs(mesh->cells[i].yp));
    }
    
    mesh->cells[i].rho=1;
    mesh->cells[i].mu=1;
    mesh->cells[i].lam=1;
    mesh->cells[i].cp=1;
    mesh->cells[i].nu=1;
    
    mesh->cells[i].k=1.0/sqrt(0.09);
    mesh->cells[i].epsilon=Cmu*pow(mesh->cells[i].k,3.0/2.0)/0.038;
    
    mesh->cells[i].Zeta=1.0/3.0* mesh->cells[i].k;
    mesh->cells[i].V2=1.0/3.0* mesh->cells[i].k;
    
    mesh->cells[i].f=1e-10;

    
    mesh->cells[i].PrT=PrT0;

    
    mesh->Senergy[i]=0.0;

    mesh->Tref = 1.0;

    mesh->rhoI_ref    =1.0;
    mesh->lambdaI_ref =1.0;
    mesh->cpI_ref     =1.0;
    mesh->muI_ref     =1.0;
    mesh->A_c         =1.0;
    mesh->B_c         =1.0;
    mesh->AsBs        =1.0;
    mesh->BsI         =1.0;

    
  }
}

/*
 * Prediction of axial velocity
 */
void Predict_w(Mesh *mesh) {
  int i;
  scalar TauFp,TauFm;
  scalar muFp,muFm;

  field l,m,r,rhs;

  scalar boundL, boundR;

  l   = mesh->l;
  m   = mesh->m;
  r   = mesh->r;
  rhs = mesh->rhs;  

  // boundary conditions
  boundL = 1.0;                                            // Neumann @ center
  boundR = - mesh->cells[imax+1].dx/mesh->cells[imax].dx;  // Dirichlet @ wall

#if VERBOSE
  printf("%f %f\n",boundL,boundR);
#endif

  for(i=1;i<=imax;i+=1) {
    //printf("%d\n",i);
    muFp=( (mesh->cells[i+1].mu+mesh->cells[i  ].mu)*0.5 + (mesh->cells[i+1].muT+mesh->cells[i  ].muT)*0.5 )/ ReTau;
    muFm=( (mesh->cells[i  ].mu+mesh->cells[i-1].mu)*0.5 + (mesh->cells[i  ].muT+mesh->cells[i-1].muT)*0.5 )/ ReTau;

    TauFp=muFp*(mesh->cells[i+1].w-mesh->cells[i  ].w)/(mesh->cells[i+1].xc-mesh->cells[i  ].xc);
    TauFm=muFm*(mesh->cells[i  ].w-mesh->cells[i-1].w)/(mesh->cells[i  ].xc-mesh->cells[i-1].xc);

    l[i-1]= -(1.0-CN) * mesh->cells[i].xFm * muFm/(mesh->cells[i].xc-mesh->cells[i-1].xc);
    
    m[i-1]= ddt*mesh->cells[i].rho +		\
      (1.0-CN) *						\
      (mesh->cells[i].xFp * muFp/(mesh->cells[i+1].xc - mesh->cells[i  ].xc) + \
       mesh->cells[i].xFm * muFm/(mesh->cells[i  ].xc - mesh->cells[i-1].xc));

    r[i-1]= -(1.0-CN) * mesh->cells[i].xFp *muFp/(mesh->cells[i+1].xc-mesh->cells[i].xc);

    rhs[i-1]= - mesh->Smom + CN*(mesh->cells[i].xFp*TauFp - mesh->cells[i].xFm*TauFm) + ddt*mesh->cells[i].rho*mesh->cells[i].w;

#if VERBOSE
    printf("%d : %f %f %f %f\n",i,l[i-1],m[i-1],r[i-1],rhs[i-1]);
#endif
  }
  ThomasAlg(imax-1,l,m,r,rhs,boundL,boundR);

  for(i=1;i<=imax;i+=1) {
    mesh->cells[i].wTs=rhs[i-1];
  }
}




/*
 * Correction of axial velocity
 */
scalar Correct_w(Mesh *mesh) {
  int i;

  scalar it,residuum; 
  
  residuum=0.0;
  for(i=1;i<=imax;i+=1) {
    it=mesh->cells[i].wTs;
    residuum+=fabs(mesh->cells[i].w-it);
    //printf("%d : %f %f %f\n",i,mesh->cells[i].w,it,residuum);
    mesh->cells[i].w   = it;

#if VERBOSE
    printf("%d : %f\n",i,mesh->cells[i].w);
#endif
    //printf("%d : %f\n",i,mesh->cells[i].w);
  }
  BoundFlow(mesh);
  

  
  return residuum;
}




/*
 * Prediction of enthalpy
 */
void Predict_h(Mesh *mesh) {
  int i;
  scalar QFp,QFm;
  scalar aFp,aFm;

  field l,m,r,rhs;

  scalar boundL, boundR;

  l   = mesh->l;
  m   = mesh->m;
  r   = mesh->r;
  rhs = mesh->rhs;  

  // boundary conditions
  boundL = 1.0;                                            // Neumann @ center
  boundR = - mesh->cells[imax+1].dx/mesh->cells[imax].dx;  // Dirichlet @ wall

#if VERBOSE
  printf("%f %f\n",boundL,boundR);
#endif

  SetEnergySourceTerm(mesh);

  for(i=1;i<=imax;i+=1) {
    //printf("%d\n",i);
    aFp=(mesh->cells[i+1].rho+mesh->cells[i  ].rho)*0.5*( (mesh->cells[i+1].a+mesh->cells[i  ].a)*0.5/PrW + (mesh->cells[i+1].aT+mesh->cells[i  ].aT)*0.5 )/ ReTau;
    aFm=(mesh->cells[i  ].rho+mesh->cells[i-1].rho)*0.5*( (mesh->cells[i  ].a+mesh->cells[i-1].a)*0.5/PrW + (mesh->cells[i  ].aT+mesh->cells[i-1].aT)*0.5 )/ ReTau;

    QFp=aFp*(mesh->cells[i+1].h/mesh->cells[i+1].cp-mesh->cells[i  ].h/mesh->cells[i  ].cp)/(mesh->cells[i+1].xc-mesh->cells[i  ].xc);
    QFm=aFm*(mesh->cells[i  ].h/mesh->cells[i  ].cp-mesh->cells[i-1].h/mesh->cells[i-1].cp)/(mesh->cells[i  ].xc-mesh->cells[i-1].xc);

    l[i-1]= -(1.0-CN) * mesh->cells[i].xFm * aFm/(mesh->cells[i].xc-mesh->cells[i-1].xc);
    
    m[i-1]= ddt*mesh->cells[i].rho +		\
      (1.0-CN) *						\
      (mesh->cells[i].xFp * aFp/(mesh->cells[i+1].xc - mesh->cells[i  ].xc) + \
       mesh->cells[i].xFm * aFm/(mesh->cells[i  ].xc - mesh->cells[i-1].xc));

    r[i-1]= -(1.0-CN) * mesh->cells[i].xFp *aFp/(mesh->cells[i+1].xc-mesh->cells[i].xc);

    rhs[i-1]= mesh->Senergy[i] + CN*(mesh->cells[i].xFp*QFp - mesh->cells[i].xFm*QFm) + ddt*mesh->cells[i].hTn;

#if VERBOSE
    printf("%d : %f %f %f %f\n",i,l[i-1],m[i-1],r[i-1],rhs[i-1]);
#endif
  }
  ThomasAlg(imax-1,l,m,r,rhs,boundL,boundR);

  for(i=1;i<=imax;i+=1) {
    mesh->cells[i].hTs=rhs[i-1];
  }
  BoundEnergyFlux(mesh);
}





/*
 * Correction of enthalpy
 */
scalar Correct_h(Mesh *mesh) {
  int i;

  scalar it,residuum;
  
  residuum=0.0;
  for(i=0;i<=imax+1;i+=1) {
    it=mesh->cells[i].hTs;
    residuum+=fabs(mesh->cells[i].h-it);
    //printf("%d : %f %f %f\n",i,mesh->cells[i].h,it,residuum);
    mesh->cells[i].h   = it;
    mesh->cells[i].hTn = mesh->cells[i].rho*mesh->cells[i].h;

#if VERBOSE
    printf("%d : %f\n",i,mesh->cells[i].h);
#endif
    //printf("%d : %f\n",i,mesh->cells[i].h);
  }

  
  return 0.0;
}

