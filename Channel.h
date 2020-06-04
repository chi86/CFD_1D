/*                                                                *
 *                    header file for flow situation              *
 *                    channel flow                                *
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
#ifndef CHANNEL_H
#define CHANNEL_H

void SetMeshMetric(Mesh *mesh);
void SetMomSourceTerm(Mesh *mesh);
void SetEnergySourceTerm(Mesh *mesh);

#endif //CHANNEL_H
