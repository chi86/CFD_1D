/****************************** turbulence model header file *******************************/
/***************************************** chi86 *******************************************/
/**************************************** Oct 2019 *****************************************/
/*******************************************************************************************/
#include "Modules.h"

#include "Utilities.h"
#include "Mesh.h"

#include "Input.h"

#ifndef TURBULENCE_H
#define TURBULENCE_H

void kEv2f(Mesh *mesh, scalar CmuKv2f_0,
	               scalar sigmaE_0,
	               scalar CE2_0,
	               scalar CL_0,
	               scalar Ceta_0,
	               scalar C1_0,
	               scalar C2_0);

#endif //TURBULENCE_H
