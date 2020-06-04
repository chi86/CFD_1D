/*                                                                *
 *                    c file for flow situation                   *
 *                    pipe flow                                   *
 *                                                                *
 *                    chi86                                       *
 *                    12.02.2020                                  *
 *                                                                *
 */

#include "Pipe.h"

/*
 * Analytical solution of the axial velocity for a laminar pipe flow
 */
scalar W_laminar_pipe(scalar ReTau,scalar mu,scalar r) {
  return ReTau/mu*(0.25-r*r);
}


/*
 * set the values for A_face/V_cell for either the F0 and F1
 */
void SetMeshMetric(Mesh *mesh) {
  int i;
  for(i=1;i<=imax;i+=1) {
    Set_x_Metric(&mesh->cells[i], \
		 mesh->cells[i  ].xu / (mesh->cells[i].xc*mesh->cells[i].dx), \
		 mesh->cells[i-1].xu / (mesh->cells[i].xc*mesh->cells[i].dx) \
		 );
  }

}

/*
 * Sourceterm for the axial momentum eqution
 */
void SetMomSourceTerm(Mesh *mesh) {
  mesh->Smom=-4.0;
}


/*
 * mass-flux (non-dimentional: mDot/(\pi*rho_0*w_tau*D^2)
 */
scalar compute_mDot(Mesh *mesh) {
  int i;

  mesh->mDot=0.0;
  for(i=1;i<=imax;i+=1) {
    mesh->mDot+=mesh->cells[i].rho*mesh->cells[i].w*mesh->cells[i].xc*mesh->cells[i].dx;
  }
  mesh->mDot*=2;
  
  return mesh->mDot;
}

/*
 * Sourceterm for the energy eqution
 */
void SetEnergySourceTerm(Mesh *mesh) {
  int i;
  scalar mDot;
  
  mDot=compute_mDot(mesh);

  for(i=1;i<=imax;i+=1) {
    mesh->Senergy[i]=mesh->cells[i].rho*mesh->cells[i].w/mDot;
  }
}
