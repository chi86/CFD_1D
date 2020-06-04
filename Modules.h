

#ifndef MODULES_H
#define MODULES_H

#include <float.h>

#ifdef USE_DOUBLES
/*
 * Singel precision
 */
typedef double scalar;
typedef double * field;

#define PREC_MIN DBL_MIN

#else
/*
 * Double precision
 */
typedef float scalar;
typedef float * field;

#define PREC_MIN FLT_MIN

#endif




#endif //MODULES_H
