#include "Cell.h"

void Init(Cell *cell) {
  cell->u=0.0;
  cell->v=0.0;
  cell->w=0.0;
  
  cell->t=0.0;
  cell->h=0.0;

  cell->uTn=0.0;
  cell->vTn=0.0;
  cell->wTn=0.0;

  cell->uTs=0.0;
  cell->vTs=0.0;
  cell->wTs=0.0;

  cell->lam=1.0;
  cell->cp=1.0;
  cell->mu=1.0;
  cell->nu=1.0;
  cell->rho=1.0;
  
  cell->a=1.0;
}

void Set_Geo_values(Cell *cell, scalar xu0,scalar xc0,scalar dx0,scalar y0,scalar yp0) {
  cell->xu=xu0;
  cell->xc=xc0;
  cell->dx=dx0;
  cell->y=y0;
  cell->yp=yp0;
}

void Set_x_Metric(Cell *cell,scalar xFp0,scalar xFm0) {
  cell->xFp=xFp0;
  cell->xFm=xFm0;
}

scalar ReturnXu(Cell *cell) { return cell->xc; }
scalar ReturnXc(Cell *cell) { return cell->xu; }
scalar ReturnDx(Cell *cell) { return cell->dx; }

scalar ReturnU(Cell *cell) { return cell->u; }
scalar ReturnV(Cell *cell) { return cell->v; }
scalar ReturnW(Cell *cell) { return cell->w; }

scalar ReturnT(Cell *cell) { return cell->t; }
scalar ReturnH(Cell *cell) { return cell->h; }
