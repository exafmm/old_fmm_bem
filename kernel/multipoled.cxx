#include <cmath>
#include <complex>

extern float *fac;
extern std::complex<double> *bnm,*bth;

void multipoled(double r, double th, double ph, int mp) {
  int m,n,nm;
  double xx,yy,s2,fact,pn,p,p1,p2;

  xx = cos(th);
  yy = sin(th);
  s2 = sqrt((1-xx)*(1+xx));
  fact = 1;
  pn = 1;
  for( m=0; m<mp; m++ ) {
    p = pn;
    nm = m*m+2*m;
    bnm[nm] = fac[nm]*p;
    p1 = p;
    p = xx*(2*m+1)*p;
    bth[nm] = fac[nm]*(p-(m+1)*xx*p1)/yy;
    for( n=m+1; n<mp; n++ ) {
      nm = n*n+n+m;
      bnm[nm] = fac[nm]*p;
      p2 = p1;
      p1 = p;
      p = (xx*(2*n+1)*p1-(n+m)*p2)/(n-m+1);
      bth[nm] = fac[nm]*((n-m+1)*p-(n+1)*xx*p1)/yy;
    }
    pn = -pn*fact*s2;
    fact += 2;
  }
}
