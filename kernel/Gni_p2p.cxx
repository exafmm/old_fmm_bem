#include "../include/constants.h"

extern float *xi,*yi,*zi,*vi,*gxi,*gyi,*gzi,*xj,*yj,*zj,*vj;
extern int **ndi,**ndj,*nej,*nij,(*neij)[nebm];
extern double *tfmm;

void Gni_p2p(int lbi, int lbj) {
  int ii,ij,jj,i,j;
  double pij,dxij,dyij,dzij,rij,rsij;
  tic = get_time();

  for( ii=0; ii<lbi; ii++ ) {
    for( ij=0; ij<nij[ii]; ij++ ) {
      jj = neij[ii][ij];
      for( i=ndi[0][ii]; i<=ndi[1][ii]; i++ ) {
        pij = 0;
        for( j=ndj[0][jj]; j<=ndj[1][jj]; j++ ) {
          dxij = xi[i]-xj[j];
          dyij = yi[i]-yj[j];
          dzij = zi[i]-zj[j];
          rij = dxij*dxij+dyij*dyij+dzij*dzij+eps;
          rsij = 0.25/pi*vj[j]/sqrt(rij)/rij;
          pij += (dxij*gxi[i]+dyij*gyi[i]+dzij*gzi[i])*rsij;
        }
        vi[i] += pij;
      }
    }
  }

  toc = get_time();
  tfmm[17] += toc-tic;
}
