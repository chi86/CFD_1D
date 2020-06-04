/*                                                                                *
 *                    CFD_1D                                                      *
 *                    1D CDF routine statistically steady flow                    *
 *                                                                                *
 *                    chi86                                                       *
 *                    12.10.2019                                                  *
 *                                                                                *
 */

#include <stdio.h>
#include <time.h>

#include "Modules.h"

#include "Mesh.h"
#include "CFD.h"

#include "Pipe.h"

#include "Fluid_oil.h"

#include "Input.h"


void Initialisation(scalar ReTau0,scalar PrW0) {
  

  
  /* ReTau=360.00;    // Reynolds-number */
  /* PrW=21.3418 ;   // moelcular Prandtl number =nuW/aW=muW*cpW/lambdaW */

  ReTau=ReTau0;    // Reynolds-number
  PrW=PrW0;   // moelcular Prandtl number =nuW/aW=muW*cpW/lambdaW
  

  imax=256;        // max cells in x
  CN=0;            // Crank Nicolson coefficient
  ddt=1;            // timestep
  RES_MIN=1E-4;    // convergence criteria
  
  PrT0=0.85;       // turbulent Prandtl number
  kappa=0.4;      // van Karman constant
  Ap=28;
  beta=5.17;      //integration constant
  Cmu=0.09;



  /*
   * Fluid properties
   * Oil data (Shell heattransfer oil S2)
   */
  T0       = 472.4601;
  w_tau    = 0.0419;
  T_tau    = 0.1948;
  h_tau    = 496.5635;
  //Density: rho(T) = A_rho + B_rho * T
  A_rho    = 1045.;
  B_rho    = -0.616;
  //Contactivity: lambda(T) = A_lambda + B_lambda * T
  A_lambda = 0.157;
  B_lambda = -7.328E-5;
  //Specific heat capacity: cp(T) = A_c + B_c * T
  A_c      = 0.818;
  B_c      = 3.664E-3;
  //Dynamic viscosity: mu(T) = A_mu * exp( B_mu / (T + C_mu) )
  A_mu     = 5.894E-5;
  B_mu     = 857.4;
  C_mu     = -172.2;




  /* // laminar */
  
  /* ReTau=100.00;    // Reynolds-number */

  
}


void Simulation() {
  int t;
  //int f;
  scalar residuum;
  
  Mesh *mesh=Mesh_new(imax,ReTau);
  printf("imax: %f\n",ReturnImax(mesh));
  printf("--------------------\n\n");

  CreateMesh(mesh);
  SetMeshMetric(mesh);
  
  InitializeFlow(mesh);
  SetMomSourceTerm(mesh);
  
#if VARMATPROPS
  RefFluidProps(mesh);
#endif
  
  residuum=1.0E0;
  
  t=0;

  //for(f=0; f<1; f+=1) {
  while(residuum>RES_MIN) {
#if !LAMINAR
    kEv2f(mesh);
#endif

    Predict_h(mesh);
    residuum=Correct_h(mesh);

#if VARMATPROPS
    hToT(mesh);
    SetFluidProperties(mesh);
#endif
    
    Predict_w(mesh);
    residuum+=Correct_w(mesh);
    
    printf("time: %5d  residuum : %.2e\n",t,residuum);
    
    //Plot_inc(mesh,t);
    
    t+=1;
  }

    
  Plot(mesh);


  FreeMesh(mesh);

}











int main(int argc, char *argv[]) {


  /****** commandline parameters ******/
  int counter; 
  printf("Program Name Is: %s",argv[0]); 
  if(argc==1) { 
    printf("\nNo Extra Command Line Argument Passed Other Than Program Name");
  }
  if(argc>=2) { 
    printf("\nNumber Of Arguments Passed: %d",argc); 
    printf("\n----Following Are The Command Line Arguments Passed----"); 
    for(counter=0;counter<argc;counter++) {
      printf("\nargv[%d]: %s",counter,argv[counter]);
    }
  }
  printf("\n\n");
  

  /****** start of main program ******/
  struct timespec before, after;
  long elapsed_secs;
  scalar ReTau0,PrW0;
  
  
  printf("\n1D CDF routine\n");
  printf("--------------\n\n");

  clock_gettime(CLOCK_REALTIME, &before);


  if(argc>2) {
    ReTau0=strtof(argv[1],NULL);
    PrW0  =strtof(argv[2],NULL);
  }
  else {
    ReTau0=360.00;    // Reynolds-number
    PrW0=21.3418 ;   // moelcular Prandtl number =nuW/aW=muW*cpW/lambdaW
  }
  printf("ReTau=%f, PrW=%f\n",ReTau0,PrW0);
  printf("\n");

 
  Initialisation(ReTau0,PrW0);
  Simulation();


  clock_gettime(CLOCK_REALTIME, &after);
  elapsed_secs = (after.tv_sec - before.tv_sec);
  printf("\nTime taken: %ld seconds\n",elapsed_secs);
  
  return 0;
}
