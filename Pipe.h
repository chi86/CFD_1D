/*                                                                *
 *                    header file for flow situation              *
 *                    pipe flow                                   *
 *                                                                *
 *                    chi86                                       *
 *                    12.02.2020                                  *
 *                                                                *
 */

#include "Modules.h"

#include "Utilities.h"
#include "Mesh.h"

#include "Input.h"


/*
 * pipe flow stuff
 */
#ifndef PIPE_H
#define PIPE_H

scalar W_laminar_pipe(scalar ReTau,scalar mu,scalar r);

void SetMeshMetric(Mesh *mesh);
void SetMomSourceTerm(Mesh *mesh);
void SetEnergySourceTerm(Mesh *mesh);

#endif //PIPE_H
