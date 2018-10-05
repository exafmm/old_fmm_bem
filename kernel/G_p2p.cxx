#include "../include/constants.h"

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi,*vi,*xj,*yj,*zj,*vj;
extern int **ndi,**ndj,*nej,*nij,(*neij)[nebm];
extern double *tfmm;

void G_p2p(int lbi, int lbj) {
  int ii,ij,jj,i,j;
  double pij,fxij,fyij,fzij,dxij,dyij,dzij,rij,rsij;
  tic = get_time();

  for( ii=0; ii<lbi; ii++ ) {
    for( ij=0; ij<nij[ii]; ij++ ) {
      jj = neij[ii][ij];
      for( i=ndi[0][ii]; i<=ndi[1][ii]; i++ ) {
        fxij = 0;
        fyij = 0;
        fzij = 0;
        pij = 0;
        for( j=ndj[0][jj]; j<=ndj[1][jj]; j++ ) {
          dxij = xi[i]-xj[j];
          dyij = yi[i]-yj[j];
          dzij = zi[i]-zj[j];
          rij = dxij*dxij+dyij*dyij+dzij*dzij+eps;
          rsij = 0.25/pi*vj[j]/rij/sqrt(rij);
          pij += 0.25/pi*vj[j]/sqrt(rij);
          fxij -= dxij*rsij;
          fyij -= dyij*rsij;
          fzij -= dzij*rsij;
        }
        gxi[i] += fxij;
        gyi[i] += fyij;
        gzi[i] += fzij;
        vi[i] += pij;
      }
    }
  }

  toc = get_time();
  tfmm[17] += toc-tic;
}
