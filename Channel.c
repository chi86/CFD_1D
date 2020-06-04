/*                                                                *
 *                    c file for flow situation                   *
 *                    channel flow                                *
 *                                                                *
 *                    chi86                                       *
 *                    12.02.2020                                  *
 *                                                                *
 */

#include "Channel.h"

/*
 * set the values for A_face/V_cell for either the F0 and F1
 */
void SetMeshMetric(Mesh *mesh) {
  int i;
  for(i=1;i<=imax;i+=1) {
    Set_x_Metric(&mesh->cells[i], \
		 1.0 / (mesh->cells[i].dx), \
		 1.0 / (mesh->cells[i].dx) \
		 );
  }

}

/*
 * Sourceterm for the axial momentum eqution
 */
void SetMomSourceTerm(Mesh *mesh) {
  mesh->Smom=-2.0;
}


/*
 * mass-flux (non-dimentional: mDot/(rho_0*w_tau*B*h/2)
 */
scalar compute_mDot(Mesh *mesh) {
  int i;

  mesh->mDot=0.0;
  for(i=1;i<=imax;i+=1) {
    mesh->mDot+=mesh->cells[i].rho*mesh->cells[i].w*mesh->cells[i].dx;
  }
  
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


    /*  def compute_cf(self): */
    /*     """compute cf""" */

    /*     imax=self.imax */

    /*     wb=0. */
    /*     rhob=0. */
    /*     for cell in self.cells[imax:0:-1]: */
    /*         wb=wb+cell.rho*cell.w*cell.rp*cell.dr */
    /*         rhob=rhob+cell.rho*cell.rp*cell.dr */

    /*     self.rhob=rhob*8.0 */
    /*     self.wb=8.0*wb/self.rhob */

    /*     self.cf=2.0/(self.rhob*self.wb**2) */
        

    /* def compute_Nusselt(self): */
    /*     """compute Nu""" */

    /*     imax=self.imax */

    /*     hb=0. */
    /*     thetab=0. */
    /*     cpb=0. */
    /*     for cell in self.cells[imax:0:-1]: */
    /*         hb=hb+cell.rho*cell.w*cell.chi*cell.rp*cell.dr */
    /*         thetab=thetab+cell.rho*cell.w*cell.theta*cell.rp*cell.dr */
    /*         cpb=cpb+cell.cp*cell.rp*cell.dr */

    /*     self.cpb=cpb*8.0 */
    /*     mdot=self.compute_mDot() */
    /*     self.hb=2.0*hb/mdot */
    /*     self.thetab=2.0*thetab/mdot */

    /*     self.Nu=self.ReTau*self.PrW/(self.thetab) */
