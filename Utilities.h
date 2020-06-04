#ifndef UTILITIES_H
#define UTILITIES_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <float.h>

#include "Modules.h"

#define max(a,b)             \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b;       \
})

#define min(a,b)             \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b;       \
})

scalar FuncRtbis(scalar x,scalar b);

scalar Rtbis(scalar b,scalar x1,scalar x2);

scalar Wlin(scalar yp);
scalar Wlog(scalar kappa,scalar beta,scalar yp);
  
scalar Hlin(scalar Pr,scalar yp);
scalar Hlog(scalar kappa,scalar beta,scalar Pr,scalar PrT,scalar yp);

void ThomasAlg(int imax, field l,field m,field r,field rhs,scalar boundL,scalar boundR);
void ThomasAlgWallModel(int imax, field l,field m,field r,field rhs,scalar boundL,scalar wallVal);

void ThomasAlgRelax(int imax, field l,field m,field r,field rhs,field x,scalar boundL,scalar boundR,scalar omega);
void ThomasAlgWallModelRelax(int imax, field l,field m,field r,field rhs,field x,scalar boundL,scalar wallVal,scalar omega);


#endif /* UTILITIES_H */
