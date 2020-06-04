#include "Mesh.h"

/*
 * Create ned mesh & allocate all arrays
 */
Mesh *Mesh_new(int imax0, scalar ReTau0) { 
  Mesh *mesh = malloc(sizeof(Mesh));

  int i;

  mesh->imax = imax0;
  mesh->ReTau = ReTau0;

  mesh->cells=malloc((imax0+2) * sizeof(Cell));

  for(i=0;i<=imax0+1;i+=1) {
    Init(&mesh->cells[i]);
  }
  
  mesh->Senergy=malloc((imax0+2) * sizeof(scalar));

  // for linear equation system solving
  mesh->l  =malloc((imax)*sizeof(scalar));
  mesh->m  =malloc((imax)*sizeof(scalar));
  mesh->r  =malloc((imax)*sizeof(scalar));
  mesh->rhs=malloc((imax)*sizeof(scalar));
  mesh->x=malloc((imax)*sizeof(scalar));
  
  /* mesh->l_k  =malloc((imax)*sizeof(scalar)); */
  /* mesh->m_k  =malloc((imax)*sizeof(scalar)); */
  /* mesh->r_k  =malloc((imax)*sizeof(scalar)); */
  /* mesh->rhs_k=malloc((imax)*sizeof(scalar)); */
  
  /* mesh->l_eps  =malloc((imax)*sizeof(scalar)); */
  /* mesh->m_eps  =malloc((imax)*sizeof(scalar)); */
  /* mesh->r_eps  =malloc((imax)*sizeof(scalar)); */
  /* mesh->rhs_eps=malloc((imax)*sizeof(scalar)); */
  
  /* mesh->l_v2  =malloc((imax)*sizeof(scalar)); */
  /* mesh->m_v2  =malloc((imax)*sizeof(scalar)); */
  /* mesh->r_v2  =malloc((imax)*sizeof(scalar)); */
  /* mesh->rhs_v2=malloc((imax)*sizeof(scalar)); */
  
  /* mesh->l_f  =malloc((imax)*sizeof(scalar)); */
  /* mesh->m_f  =malloc((imax)*sizeof(scalar)); */
  /* mesh->r_f  =malloc((imax)*sizeof(scalar)); */
  /* mesh->rhs_f=malloc((imax)*sizeof(scalar)); */
  
  return mesh;
}

/*
 * Erase all allocated arrays
 */
void FreeMesh(Mesh *mesh) {
  free(mesh->cells);
  
  free(mesh->l  );
  free(mesh->m  );
  free(mesh->r  );
  free(mesh->rhs);
  free(mesh->x);
  
  /* free(mesh->l_k  ); */
  /* free(mesh->m_k  ); */
  /* free(mesh->r_k  ); */
  /* free(mesh->rhs_k); */
  
  /* free(mesh->l_eps  ); */
  /* free(mesh->m_eps  ); */
  /* free(mesh->r_eps  ); */
  /* free(mesh->rhs_eps); */
  
  /* free(mesh->l_v2  ); */
  /* free(mesh->m_v2  ); */
  /* free(mesh->r_v2  ); */
  /* free(mesh->rhs_v2); */
  
  /* free(mesh->l_f  ); */
  /* free(mesh->m_f  ); */
  /* free(mesh->r_f  ); */
  /* free(mesh->rhs_f); */
  
  free(mesh->Senergy);
  
  free(mesh);
}
  

/*
 * Mesh member functions
 */
scalar ReturnImax(Mesh *mesh) { return mesh->imax; }

/*
 * Create the computational mesh
 */
void CreateMesh(Mesh *mesh) {
  scalar DetaW,DetaC;
  scalar A,B;
  scalar ppt;

  scalar *Ru;
  scalar *Rp;
  scalar *dr;
  
  scalar eta,dr_rel;
  
  int i;

  Ru=malloc((imax+2)*sizeof(scalar));
  Rp=malloc((imax+2)*sizeof(scalar));
  dr=malloc((imax+2)*sizeof(scalar));

  // Vinokur, Marcel (1983)
  // On one-dimensional stretching functions for finite-difference calculations.
  // Journal of Computational Physics, 50, 215-234.

  // Parameters
  // DetaW ... Delta y^+ @ wall
  DetaW = 0.05;
  // DetaC ... Delta y^+ @ centre
  DetaC = 0.40;
    
  // computation
  DetaW = DetaW*2.0/mesh->ReTau;
  DetaC = DetaC*2.0/mesh->ReTau;

  A=sqrt(DetaW)/sqrt(DetaC);
  B=1.0/(imax*sqrt(DetaC*DetaW));

  ppt=Rtbis(B,0.1,20.0);

#if VERBOSE
  printf("B: %f\n",B);
  printf("ppt : %f\n",ppt);
#endif
  
  for (i=0;i<imax;i+=1) {
      eta=(scalar)i/(scalar)imax;
      dr_rel=0.5*(1.0+tanh(ppt*(eta-0.5))/tanh(ppt/2.0));

      Ru[i]=0.5*dr_rel/(A+(1.0-A)*dr_rel);
    }

  Ru[0]    = 0.0;
  Ru[imax] = 0.5;

  for (i=1;i<=imax;i+=1) {
      Rp[i] = ( Ru[i] + Ru[i-1] ) / 2.0;
      dr[i] = ( Ru[i] - Ru[i-1] );
    }
  dr[imax+1] = dr[imax];
  Ru[imax+1] = Ru[imax] + dr[imax+1];
  Rp[imax+1] = Ru[imax] + dr[imax+1] / 2.0;
  dr[0]  = dr[1];
  Rp[0]  = Ru[0] - dr[0] / 2.0;

  BoundFlow(mesh);
  BoundEnergy(mesh);
  
  for(i=0;i<=imax+1;i+=1) {
    Set_Geo_values(&mesh->cells[i],Ru[i],Rp[i],dr[i],(0.5-Rp[i]),(0.5-Rp[i])*ReTau);
    }


  mesh->EWall=2.0/(powf(mesh->cells[imax].y,2.0))/ReTau;

  

  free(Ru);
  free(Rp);
  free(dr);
}


/*
 * setting flow boundary condition
 */
void BoundFlow(Mesh *mesh) {
  mesh->cells[0     ].w =   mesh->cells[1].w;
  mesh->cells[imax+1].w = - mesh->cells[imax].w;
}

/*
 * setting flow boundary condition
 */
void BoundEnergy(Mesh *mesh) {
  mesh->cells[0     ].h =   mesh->cells[1].h;
  mesh->cells[imax+1].h = - mesh->cells[imax].h;
}

/*
 * setting flow boundary condition for fluxes
 */
void BoundFlowFlux(Mesh *mesh) {
  mesh->cells[0     ].wTs =   mesh->cells[1].wTs;
  mesh->cells[imax+1].wTs = - mesh->cells[imax].wTs;
}

/*
 * setting flow boundary condition for fluxes
 */
void BoundEnergyFlux(Mesh *mesh) {
  mesh->cells[0     ].hTs =   mesh->cells[1].hTs;
  mesh->cells[imax+1].hTs = - mesh->cells[imax].hTs;
}


// for the pipe case
void Plot(Mesh *mesh) {
  int i;

  FILE *write_ptr;

  ComputeTau(mesh);
  ComputeQ(mesh);
  

  write_ptr = fopen("data.dat","w");

  fprintf(write_ptr,"%20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s\n","ru","rp","dr","yp","W","T","H","k","epsilon","Zeta","f","muT","L","T","aT","tau_l","tau_t","q_l","q_t","V2");
  for(i=1;i<=imax;i=i+1)
    {
      fprintf(write_ptr,"%20.12e %20.12e %20.12e %20.12e ",mesh->cells[i].xu,mesh->cells[i].xc,mesh->cells[i].dx,mesh->cells[i].yp);
      fprintf(write_ptr,"%20.12e ",mesh->cells[i].w);
      fprintf(write_ptr,"%20.12e %20.12e",mesh->cells[i].t,mesh->cells[i].h);
      fprintf(write_ptr,"%20.12e %20.12e %20.12e %20.12e %20.12e",mesh->cells[i].k,mesh->cells[i].epsilon,mesh->cells[i].Zeta,mesh->cells[i].f,mesh->cells[i].muT);
      fprintf(write_ptr,"%20.12e %20.12e",mesh->cells[i].L,mesh->cells[i].T);
      fprintf(write_ptr,"%20.12e",mesh->cells[i].aT);
      fprintf(write_ptr,"%20.12e %20.12e",mesh->cells[i].tau_l,mesh->cells[i].tau_t);
      fprintf(write_ptr,"%20.12e %20.12e",mesh->cells[i].q_l,mesh->cells[i].q_t);
      fprintf(write_ptr,"%20.12e",mesh->cells[i].V2);
      fprintf(write_ptr,"\n");
    }
  
  fclose(write_ptr);
  
}

void ComputeTau(Mesh *mesh) {
  int i;
  
  scalar muL,muT;
  scalar tauFm,tauFp;
  
  for(i=1;i<=imax;i+=1) {
    muL=(mesh->cells[i  ].mu +mesh->cells[i-1].mu)*0.5;
    tauFm=muL*(mesh->cells[i  ].w-mesh->cells[i-1].w)/(mesh->cells[i  ].xc-mesh->cells[i-1].xc);
    
    muL=(mesh->cells[i+1].mu +mesh->cells[i  ].mu)*0.5;
    tauFp=muL*(mesh->cells[i+1].w-mesh->cells[i  ].w)/(mesh->cells[i+1].xc-mesh->cells[i  ].xc);

    mesh->cells[i].tau_l=0.5*(tauFm+tauFp)/ReTau;
    
    muT=(mesh->cells[i  ].muT+mesh->cells[i-1].muT)*0.5;
    tauFm=muT*(mesh->cells[i  ].w-mesh->cells[i-1].w)/(mesh->cells[i  ].xc-mesh->cells[i-1].xc);
    
    muT=(mesh->cells[i+1].muT+mesh->cells[i  ].muT)*0.5;
    tauFp=muT*(mesh->cells[i+1].w-mesh->cells[i  ].w)/(mesh->cells[i+1].xc-mesh->cells[i  ].xc);
        
    mesh->cells[i].tau_t=0.5*(tauFm+tauFp)/ReTau;
  }
}


void ComputeQ(Mesh *mesh) {
  int i;
  
  scalar aL,aT;
  scalar qFm,qFp;
  
  for(i=1;i<=imax;i+=1) {
    aL=(mesh->cells[i  ].rho+mesh->cells[i-1].rho)*0.5*(mesh->cells[i  ].a +mesh->cells[i-1].a)*0.5;
    qFm=aL*(mesh->cells[i  ].h/mesh->cells[i  ].cp-mesh->cells[i-1].h/mesh->cells[i-1].cp)/(mesh->cells[i  ].xc-mesh->cells[i-1].xc);
    
    aL=(mesh->cells[i+1].rho+mesh->cells[i  ].rho)*0.5*(mesh->cells[i+1].a +mesh->cells[i  ].a)*0.5;
    qFp=aL*(mesh->cells[i+1].h/mesh->cells[i+1].cp-mesh->cells[i  ].h/mesh->cells[i  ].cp)/(mesh->cells[i+1].xc-mesh->cells[i  ].xc);

    mesh->cells[i].q_l=0.5*(qFm+qFp)/(ReTau*PrW);

        
    aT=(mesh->cells[i  ].rho+mesh->cells[i-1].rho)*0.5*(mesh->cells[i  ].aT+mesh->cells[i-1].aT)*0.5;
    qFm=aT*(mesh->cells[i  ].h/mesh->cells[i  ].cp-mesh->cells[i-1].h/mesh->cells[i-1].cp)/(mesh->cells[i  ].xc-mesh->cells[i-1].xc);
    
    aT=(mesh->cells[i+1].rho+mesh->cells[i  ].rho)*0.5*(mesh->cells[i+1].aT+mesh->cells[i  ].aT)*0.5;
    qFp=aT*(mesh->cells[i+1].h/mesh->cells[i+1].cp-mesh->cells[i  ].h/mesh->cells[i  ].cp)/(mesh->cells[i+1].xc-mesh->cells[i  ].xc);

    mesh->cells[i].q_t=0.5*(qFm+qFp)/(ReTau);
  }
}
