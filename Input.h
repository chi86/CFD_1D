/*                                                                *
 *                    Input data for CFD simulation               *
 *                    parameter definition                        *
 *                                                                *
 *                    chi86                                       *
 *                    12.02.2020                                  *
 *                                                                *
 */




#include "Modules.h"

#ifndef INPUT_H
#define INPUT_H

// max cells in x
int imax;

// Crank Nicolson coefficient
scalar CN;

// timestep
scalar ddt;

// convergence criteria
scalar RES_MIN;

// Reynolds-number
scalar ReTau;

// moelcular Prandtl number =nuW/aW=muW*cpW/lambdaW
scalar PrW;

// turbulent Prandtl number
scalar PrT0;

// van Karman constant
scalar kappa;

scalar Ap;

//integration constant
scalar beta;

scalar Cmu;


// fluid properties
scalar T0;
scalar w_tau;
scalar T_tau;
scalar h_tau;

scalar A_rho;
scalar B_rho;

scalar A_lambda;
scalar B_lambda;

scalar A_c;
scalar B_c;

scalar A_mu;
scalar B_mu;
scalar C_mu;




#endif /* INPUT_H */
