#include "../include/constants.h"

extern float *fac;
extern std::complex<double> *bnm,*bth;

extern void in(double, double, const int, double*, int);

void multipoleh(double rh, double al, double be, int mp) {
  int m,n,nm=0;
  double xx,s2,fact,pn,p,p1,p2,inr[mpcmp];

  in(scale,rh*lambda,mp,inr,nm);
  xx = cos(al);
  s2 = sqrt((1-xx)*(1+xx));
  fact = 1;
  pn = 1;
  for( m=0; m<mp; m++ ) {
    p = pn;
    nm = m*m+2*m;
    bnm[nm] = inr[m]*(2*m+1)*fac[nm]*fac[nm]*p;
    p1 = p;
    p = xx*(2*m+1)*p;
    for( n=m+1; n<mp; n++ ) {
      nm = n*n+n+m;
      bnm[nm] = inr[n]*(2*n+1)*fac[nm]*fac[nm]*p;
      p2 = p1;
      p1 = p;
      p = (xx*(2*n+1)*p1-(n+m)*p2)/(n-m+1);
    }
    pn = -pn*fact*s2;
    fact += 2;
  }
}
