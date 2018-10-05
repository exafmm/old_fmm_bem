#include "../include/constants.h"

extern float *xi,*yi,*zi,*vi;
extern float *xj,*yj,*zj,*gxj,*gyj,*gzj,*vj;

void Gnj_dir(int n0, int n1, int n2, int n3) {
  int i,j;
  double pij,dxij,dyij,dzij,rij,rsij;
  for( i=n0; i<=n1; i++ ) {
    pij = 0;
    for( j=n2; j<=n3; j++ ){
      dxij = xi[i]-xj[j];
      dyij = yi[i]-yj[j];
      dzij = zi[i]-zj[j];
      rij = dxij*dxij+dyij*dyij+dzij*dzij+eps;
      rsij = 0.25/pi*vj[j]/rij/sqrt(rij);
      pij += (dxij*gxj[j]+dyij*gyj[j]+dzij*gzj[j])*rsij;
    }
    vi[i] = pij;
  }
}
