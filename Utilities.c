// Utility funcitons
#include "Utilities.h"

scalar FuncRtbis(scalar x,scalar b)
{
  return sinh(x)/x-b;
}

scalar Rtbis(scalar b,scalar x1,scalar x2)
{
  int jmax=100;
  scalar rtbis;
  scalar xacc=0e-8;

  int j;
  scalar dx,f,fmid,xmid;

  fmid=FuncRtbis(x2,b);
  f=FuncRtbis(x1,b);
  if(f*fmid>=0.0)
    {
      printf("root must be bracketed in rtbis");
      return 0.0;
    }
  if(f<0.)
    {
      rtbis=x1;
      dx=x2-x1;
    }
  else
    {
      rtbis=x2;
      dx=x1-x2;
    }
  for(j=1;j<=jmax;j+=1)
    {
      dx=dx*0.5;
      xmid=rtbis+dx;
      fmid=FuncRtbis(xmid,b);
      if(fmid<=0.0)
	rtbis=xmid;
      if(abs(dx)<xacc && fmid==0.0)
	return rtbis;
    }
  return rtbis;
}

/* Axial velocity linear law */
scalar Wlin(scalar yp) {
  return yp;
}

/* Axial velocity low law */
scalar Wlog(scalar kappa,scalar beta,scalar yp) {
  return 1/kappa*log(yp)+beta;
}

/* Axial velocity linear law */
scalar Hlin(scalar Pr,scalar yp) {
  return Pr*yp;
}

/* Axial velocity low law */
scalar Hlog(scalar kappa,scalar beta,scalar Pr,scalar PrT,scalar yp) {
  scalar P;
  P=9.24*( pow(Pr/PrT,0.75)-1 )*(1+0.28*exp(-0.007*(Pr/PrT)));
  return PrT*(1/kappa*log(yp)+beta+P);
}

/* Thomas algorithm for solving a linear equation system of the form: A.w=rhs, A=l*w-1+m*w0+r*w+1 */
void ThomasAlg(int imax, field l,field m,field r,field rhs,scalar boundL,scalar boundR) {

  int i;

  // set boundary conditions
  m[0]    = m[0]   +l[0]   *boundL;
  m[imax] = m[imax]+r[imax]*boundR;

  // i=0
  r[0]=r[0]/m[0];
  rhs[0]=rhs[0]/m[0];

  //i=1:imax-1
  for(i=1;i<=imax-1;i+=1) {
    r[i] = r[i]/(m[i]-r[i-1]*l[i]);
    rhs[i] = (rhs[i]-rhs[i-1]*l[i])/(m[i]-r[i-1]*l[i]);
  }

  // i=imax;
  rhs[imax] = (rhs[imax]-rhs[imax-1]*l[imax])/(m[imax]-r[imax-1]*l[imax]);

  // i=imax-1:0
  for(i=imax-1;i>=0;i-=1) {
    //printf("%d : %f\n",i,rhs[i]);
    rhs[i] = rhs[i] - r[i] * rhs[i+1];
  }
}

/* Thomas algorithm for solving a linear equation system of the form: A.w=rhs, A=l*w-1+m*w0+r*w+1 */
void ThomasAlgWallModel(int imax, field l,field m,field r,field rhs,scalar boundL,scalar wallVal) {

  int i;

  // set boundary conditions
  m[0]    = m[0]   +l[0]   *boundL;

  // i=0
  r[0]=r[0]/m[0];
  rhs[0]=rhs[0]/m[0];

  //i=1:imax-1
  for(i=1;i<=imax-1;i+=1) {
    r[i] = r[i]/(m[i]-r[i-1]*l[i]);
    rhs[i] = (rhs[i]-rhs[i-1]*l[i])/(m[i]-r[i-1]*l[i]);
  }

  // i=imax;
  rhs[imax] = wallVal;

  // i=imax-1:0
  for(i=imax-1;i>=0;i-=1) {
    //printf("%d : %f\n",i,rhs[i]);
    rhs[i] = rhs[i] - r[i] * rhs[i+1];
  }
}



/* Thomas algorithm for solving a linear equation system of the form: A.w=rhs, A=l*w-1+m*w0+r*w+1 */
void ThomasAlgRelax(int imax, field l,field m,field r,field rhs,field x,scalar boundL,scalar boundR,scalar omega) {

  int i;

  // set boundary conditions
  m[0]    = m[0]   +l[0]   *boundL;
  m[imax] = m[imax]+r[imax]*boundR;

  //under-relaxation
  for(i=0;i<=imax;i+=1) {
    rhs[i]=rhs[i] + (1.0-omega)/omega*m[i]*x[i];
    m[i]=m[i]/omega;
  }

  // i=0
  r[0]=r[0]/m[0];
  rhs[0]=rhs[0]/m[0];

  //i=1:imax-1
  for(i=1;i<=imax-1;i+=1) {
    r[i] = r[i]/(m[i]-r[i-1]*l[i]);
    rhs[i] = (rhs[i]-rhs[i-1]*l[i])/(m[i]-r[i-1]*l[i]);
  }

  // i=imax;
  rhs[imax] = (rhs[imax]-rhs[imax-1]*l[imax])/(m[imax]-r[imax-1]*l[imax]);

  // i=imax-1:0
  for(i=imax-1;i>=0;i-=1) {
    //printf("%d : %f\n",i,rhs[i]);
    rhs[i] = rhs[i] - r[i] * rhs[i+1];
  }
}


/* Thomas algorithm for solving a linear equation system of the form: A.w=rhs, A=l*w-1+m*w0+r*w+1 */
void ThomasAlgWallModelRelax(int imax, field l,field m,field r,field rhs,field x,scalar boundL,scalar wallVal,scalar omega) {

  int i;

  // set boundary conditions
  m[0]    = m[0]   +l[0]   *boundL;

  //under-relaxation
  for(i=0;i<=imax;i+=1) {
    rhs[i]=rhs[i] + (1.0-omega)/omega*m[i]*x[i];
    m[i]=m[i]/omega;
  }

  // i=0
  r[0]=r[0]/m[0];
  rhs[0]=rhs[0]/m[0];

  //i=1:imax-1
  for(i=1;i<=imax-1;i+=1) {
    r[i] = r[i]/(m[i]-r[i-1]*l[i]);
    rhs[i] = (rhs[i]-rhs[i-1]*l[i])/(m[i]-r[i-1]*l[i]);
  }

  // i=imax;
  rhs[imax] = wallVal;

  // i=imax-1:0
  for(i=imax-1;i>=0;i-=1) {
    //printf("%d : %f\n",i,rhs[i]);
    rhs[i] = rhs[i] - r[i] * rhs[i+1];
  }
}
