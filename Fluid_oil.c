#include "Fluid_oil.h"

/*
 * reference values for the fluidproperties
 */
void RefFluidProps(Mesh *mesh) {
  scalar As_c,Bs_c;

  mesh->Tref=T0/T_tau;
  
  //density in [kg/m3]
  mesh->rhoI_ref =    1.0 / (A_rho +B_rho *T0 );
  //conductivity in [W/mK]
  mesh->lambdaI_ref = 1.0 / (A_lambda +B_lambda *T0 );
  //integral specific heat capacity in [J/kgK]
  mesh->cpI_ref =     1.0 / ( (A_c +B_c *T0) * 1000.0 );
  //dynamic viscosity in [Pas]
  mesh->muI_ref =     1.0 / (A_mu * exp(B_mu / (T0 +C_mu) ) );
  
  //non-dimensional coefficients for cp* in [J/kgK]
  As_c =A_c * 1000.0 * mesh->cpI_ref;
  Bs_c =B_c * 1000.0 * mesh->cpI_ref * T_tau;
  
  
  mesh->AsBs = As_c / Bs_c;
  mesh->BsI  = 1.0  /  Bs_c ;

}


/*
 * convert the enthalpy (in this case the non-dimensional) to temperature
 */
void hToT(Mesh *mesh) {
  int i;

  for(i=0;i<=imax+1;i+=1) {
    mesh->cells[i].t= ( -mesh->AsBs +		       \
			sqrt( pow(mesh->Tref + mesh->AsBs,2) -	\
			      2.0*mesh->cells[i].h * mesh->BsI	\
			      )					\
			)* T_tau;
    //printf("%d %f %f %f %f\n",i,mesh->cells[i].t,mesh->cells[i].h,mesh->AsBs, sqrt( pow(mesh->Tref + mesh->AsBs,2)-2.0*mesh->cells[i].h * mesh->BsI));
  }
  
}

/*
 * Compute temperature dependent thermophysical material properties of the fluid
 */
void SetFluidProperties(Mesh *mesh) {
  int i;

  for(i=0;i<=imax+1;i+=1) {
    mesh->cells[i].rho = (A_rho + B_rho * mesh->cells[i].t )* mesh->rhoI_ref;
    
    mesh->cells[i].lam = (A_lambda + B_lambda * mesh->cells[i].t ) * mesh->lambdaI_ref;

    mesh->cells[i].cp  = (A_c + B_c * mesh->cells[i].t ) * 1000. * mesh->cpI_ref;
            
    mesh->cells[i].mu  =  (A_mu * exp(B_mu/(mesh->cells[i].t+C_mu)) )* mesh->muI_ref;

    mesh->cells[i].nu = mesh->cells[i].mu/mesh->cells[i].rho;
    mesh->cells[i].a  = mesh->cells[i].lam / (mesh->cells[i].rho * mesh->cells[i].cp);
  }
  
}
