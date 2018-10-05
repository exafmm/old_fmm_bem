#include "../include/constants.h"

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi,*vi;
extern int **ndi,*nfi;
extern std::complex<double> (*ax)[mpsym],*bnm,*bth;
extern double *tfmm;

extern void boxc(int, int, int*);
extern void multipoled(double, double, double, int);

void Gn_l2p(int nmp, int mp, int lbi) {
  int kl,ii,i,n,nm,nms,m,nc[3];
  double rb,xic,yic,zic,xiic,yiic,ziic,r,th,ph,gx,gy,gz,gr,gth,gph;
  std::complex<double> aax[mpsym],rr,rth,rph,cnm(0.0,1.0),eim;
  tic = get_time();

  kl = 1 << lmax;
  rb = rd/kl;
  for( ii=0; ii<lbi; ii++ ) {
    boxc(nfi[ii],3,nc);
    xic = xmin+(nc[0]+0.5)*rb;
    yic = ymin+(nc[1]+0.5)*rb;
    zic = zmin+(nc[2]+0.5)*rb;
    for( i=0; i<nmp; i++ ) aax[i] = ax[ii][i];
    for( i=ndi[0][ii]; i<=ndi[1][ii]; i++ ) {
      xiic = xi[i]-xic;
      yiic = yi[i]-yic;
      ziic = zi[i]-zic;
      r = sqrt(xiic*xiic+yiic*yiic+ziic*ziic)+eps;
      th = acos(ziic/r);
      if( std::abs(xiic)+std::abs(yiic) < eps ) {
        ph = 0;
      } else if( std::abs(xiic) < eps ) {
        ph = yiic/std::abs(yiic)*pi*0.5;
      } else if( xiic > 0 ) {
        ph = atan(yiic/xiic);
      } else {
        ph = atan(yiic/xiic)+pi;
      }
      gr = 0;
      gth = 0;
      gph = 0;
      multipoled(r,th,ph,mp);
      for( n=0; n<mp; n++ ) {
        nm = n*n+n;
        nms = n*(n+1)/2;
        rr = n*pow(r,n-1)*bnm[nm];
        rth = pow(r,n)*bth[nm];
        gr += real(rr*aax[nms]);
        gth += real(rth*aax[nms]);
        for( m=1; m<=n; m++ ) {
          nm = n*n+n+m;
          nms = n*(n+1)/2+m;
          eim = exp(m*ph*cnm);
          rr = n*pow(r,n-1)*bnm[nm]*eim;
          rth = pow(r,n)*bth[nm]*eim;
          rph = m*pow(r,n)*bnm[nm]*eim*cnm;
          gr += 2*real(rr*aax[nms]);
          gth += 2*real(rth*aax[nms]);
          gph += 2*real(rph*aax[nms]);
        }
      }
      gx = sin(th)*cos(ph)*gr+cos(th)*cos(ph)/r*gth-sin(ph)/r/sin(th)*gph;
      gy = sin(th)*sin(ph)*gr+cos(th)*sin(ph)/r*gth+cos(ph)/r/sin(th)*gph;
      gz = cos(th)*gr-sin(th)/r*gth;
      vi[i] -= 0.25/pi*(gxi[i]*gx+gyi[i]*gy+gzi[i]*gz);
    }
  }

  toc = get_time();
  tfmm[22] += toc-tic;
}
