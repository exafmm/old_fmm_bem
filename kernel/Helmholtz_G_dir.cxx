#include "../include/constants.h"

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi,*vi;
extern float *xj,*yj,*zj,*vj;

void Helmholtz_G_dir(int n0, int n1, int n2, int n3) {
  int i,j;
  double pij,dxij,dyij,dzij,rij;
  for( i=n0; i<=n1; i++ ) {
    pij = 0;
    for( j=n2; j<=n3; j++ ){
      dxij = xi[i]-xj[j];
      dyij = yi[i]-yj[j];
      dzij = zi[i]-zj[j];
      rij = sqrt(dxij*dxij+dyij*dyij+dzij*dzij+eps);
      pij += 0.5*pi*vj[j]/rij/lambda*exp(-lambda*rij);
    }
    vi[i] = pij;
  }
}
