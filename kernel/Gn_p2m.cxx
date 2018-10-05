#include "../include/constants.h"

extern float *xj,*yj,*zj,*gxj,*gyj,*gzj,*vj;
extern int **ndj,*nfj;
extern std::complex<double> (*bx)[mpsym],*bnm,*bth;
extern double *tfmm;

extern void boxc(int, int, int*);
extern void multipoled(double, double, double, int);

void Gn_p2m(int nmp, int mp, int lbj) {
  int kl,jj,j,n,m,nm,nms,nc[3];
  double rb,xjc,yjc,zjc,xjjc,yjjc,zjjc,rh,al,be;
  std::complex<double> bbx[mpsym],brh,bal,bbe,bxd,byd,bzd,cnm(0.0,1.0),eim;
  tic = get_time();

  kl = 1 << lmax;
  rb = rd/kl;
  for( jj=0; jj<lbj; jj++ ) {
    boxc(nfj[jj],3,nc);
    xjc = xmin+(nc[0]+0.5)*rb;
    yjc = ymin+(nc[1]+0.5)*rb;
    zjc = zmin+(nc[2]+0.5)*rb;
    for( j=0; j<nmp; j++ ) {
      bbx[j] = 0;
    }
    for( j=ndj[0][jj]; j<=ndj[1][jj]; j++ ) {
      xjjc = xj[j]-xjc;
      yjjc = yj[j]-yjc;
      zjjc = zj[j]-zjc;
      rh = sqrt(xjjc*xjjc+yjjc*yjjc+zjjc*zjjc)+eps;
      al = acos(zjjc/rh);
      if( std::abs(xjjc)+std::abs(yjjc) < eps ) {
        be = 0;
      } else if( std::abs(xjjc) < eps ) {
        be = yjjc/std::abs(yjjc)*pi*0.5;
      } else if( xjjc > 0 ) {
        be = atan(yjjc/xjjc);
      } else {
        be = atan(yjjc/xjjc)+pi;
      }
      multipoled(rh,al,-be,mp);
      for( n=0; n<mp; n++ ) {
        for( m=0; m<=n; m++ ) {
          nm = n*n+n+m;
          nms = n*(n+1)/2+m;
          eim = exp(-m*be*cnm);
          brh = n*pow(rh,n-1)*bnm[nm]*eim;
          bal = pow(rh,n)*bth[nm]*eim;
          bbe = -m*pow(rh,n)*bnm[nm]*eim*cnm;
          bxd = sin(al)*cos(be)*brh+cos(al)*cos(be)/rh*bal-sin(be)/rh/sin(al)*bbe;
          byd = sin(al)*sin(be)*brh+cos(al)*sin(be)/rh*bal+cos(be)/rh/sin(al)*bbe;
          bzd = cos(al)*brh-sin(al)/rh*bal;
          bbx[nms] += ((std::complex<double>) (vj[j]*gxj[j]))*bxd;
          bbx[nms] += ((std::complex<double>) (vj[j]*gyj[j]))*byd;
          bbx[nms] += ((std::complex<double>) (vj[j]*gzj[j]))*bzd;
        }
      }
    }
    for( j=0; j<nmp; j++ ) {
      bx[jj][j] = bbx[j];
    }
  }

  toc = get_time();
  tfmm[18] += toc-tic;
}
