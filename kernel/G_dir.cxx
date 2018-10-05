#include "../include/constants.h"

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi,*vi;
extern float *xj,*yj,*zj,*vj;

void G_dir(int n0, int n1, int n2, int n3) {
  int i,j;
  double pij,fxij,fyij,fzij,dxij,dyij,dzij,rij,rsij;
  for( i=n0; i<=n1; i++ ) {
    pij = 0;
    fxij = 0;
    fyij = 0;
    fzij = 0;
    for( j=n2; j<=n3; j++ ){
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
    gxi[i] = fxij;
    gyi[i] = fyij;
    gzi[i] = fzij;
    vi[i] = pij;
  }
}
