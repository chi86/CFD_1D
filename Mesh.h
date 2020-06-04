/****************************** Mesh header file *******************************/
/*********************************** chi86 *************************************/
/********************************** Oct 2019 ***********************************/
/*******************************************************************************/
#include "Modules.h"

#include "Utilities.h"
#include "Cell.h"

#include "Input.h"

#ifndef MESH_H
#define MESH_H

struct Mesh {
  int imax;
  Cell *cells;
  
  field l,m,r,rhs,x;
  
  /* field l_k,m_k,r_k,rhs_k; */
  /* field l_eps,m_eps,r_eps,rhs_eps; */
  /* field l_v2,m_v2,r_v2,rhs_v2; */
  /* field l_f,m_f,r_f,rhs_f; */

  scalar ReTau;

  scalar Smom;
  field Senergy;

  scalar mDot;

  scalar EWall;

  /*
   *  Fluidproperties
   */
  scalar Tref;
  
  scalar rhoI_ref;
  scalar lambdaI_ref;
  scalar cpI_ref;
  scalar muI_ref;
  scalar A_c;
  scalar B_c;
  scalar AsBs;
  scalar BsI;
  
};
typedef struct Mesh Mesh;

// Constructor
Mesh *Mesh_new(int imax0, scalar ReTau0);
// Destructor
void FreeMesh(Mesh *mesh);

scalar ReturnImax(Mesh *mesh);
void CreateMesh(Mesh *mesh);
void BoundFlow(Mesh *mesh);
void BoundEnergy(Mesh *mesh);
void BoundFlowFlux(Mesh *mesh);
void BoundEnergyFlux(Mesh *mesh);

void ComputeTau(Mesh *mesh);
void ComputeQ(Mesh *mesh);

void InitializeFlow(Mesh *mesh);

void Plot(Mesh *mesh);

// public:
//   Mesh(int,scalar); // constructor
//   ~Mesh();           // destructor 

#endif //MESH_H
