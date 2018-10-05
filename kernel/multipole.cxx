#include <cmath>
#include <complex>

extern float *fac;
extern std::complex<double> *bnm,*bth;

void multipole(double rh, double al, double be, int mp) {
  int m,n,nm;
  double xx,s2,fact,pn,p,p1,p2,rhm,rhn;

  xx = cos(al);
  s2 = sqrt((1-xx)*(1+xx));
  fact = 1;
  pn = 1;
  rhm = 1;
  for( m=0; m<mp; m++ ) {
    p = pn;
    nm = m*m+2*m;
    bnm[nm] = rhm*fac[nm]*p;
    p1 = p;
    p = xx*(2*m+1)*p;
    rhm *= rh;
    rhn = rhm;
    for( n=m+1; n<mp; n++ ) {
      nm = n*n+n+m;
      bnm[nm] = rhn*fac[nm]*p;
      p2 = p1;
      p1 = p;
      p = (xx*(2*n+1)*p1-(n+m)*p2)/(n-m+1);
      rhn *=rh;
    }
    pn = -pn*fact*s2;
    fact += 2;
  }
}
