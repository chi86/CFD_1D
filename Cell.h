/****************************** Cell header file *******************************/
/*********************************** chi86 *************************************/
/********************************** Oct 2019 ***********************************/
/*******************************************************************************/
#include "Modules.h"

#ifndef CELL_H
#define CELL_H

struct Cell {
  scalar xu,xc,dx,y,yp; // geometry
  scalar u,v,w;         // velocities
  scalar uTn,vTn,wTn;   // du/dt,dv/dt,dw/dt @ n
  scalar uTs,vTs,wTs;   // du/dt,dv/dt,dw/dt @ star
  
  scalar h;             // enthalpy
  scalar hTn;   // du/dt,dv/dt,dw/dt @ n
  scalar hTs;   // du/dt,dv/dt,dw/dt @ star
  
  scalar t;             // temperature

  scalar lam;           // thermal conductivity
  scalar cp;            // specific heat capacity
  scalar mu;            // dyn. viscosity
  scalar nu;            // kin. viscosity
  scalar rho;           // density
  scalar a;             // Thermal diffusivity

  scalar strainR,PExpTurb;

  scalar tau_l,tau_t;   // shear stress
  scalar q_l,q_t;       // heat flux

  
  //metric
  scalar xFp,xFm;

  // turbulence properties
  scalar muT;
  scalar aT;
  scalar PrT;

  scalar k;
  scalar epsilon;
  scalar Zeta;
  scalar V2;
  scalar f;
  scalar CE1;
  scalar T;
  scalar L;

};
typedef struct Cell Cell;

// function prototype
void Init(Cell *cell);

void Set_Geo_values(Cell *cell, scalar ru0,scalar rp0,scalar dr0,scalar y0,scalar yp0);
void Set_x_Metric(Cell *cell,scalar xFp0,scalar xFm0);

scalar ReturnXu(Cell *cell);
scalar ReturnXc(Cell *cell);
scalar ReturnDx(Cell *cell);

scalar ReturnU(Cell *cell);
scalar ReturnV(Cell *cell);
scalar ReturnW(Cell *cell);

scalar ReturnH(Cell *cell);
scalar ReturnT(Cell *cell);


#endif //CELL_H
