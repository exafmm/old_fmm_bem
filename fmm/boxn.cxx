#include "../include/constants.h"

extern int *nfn;

void boxn(int n, float *x, float *y, float *z) {
  int kl,i,j,nf,nfx[3];
  double rb;
  kl = 1 << lmax;
  rb = kl/rd;

  for( j=0; j<n; j++ ) {
    nfx[0] = int((x[j]-xmin)*rb);
    nfx[1] = int((y[j]-ymin)*rb);
    nfx[2] = int((z[j]-zmin)*rb);
    for( i=0; i<3; i++ ) {
      if( nfx[i] >= (1 << lmax) ) nfx[i]--;
    }
    nf = 0;
    for( i=0; i<lmax; i++ ) {
      nf += nfx[0]%2 << (3*i+1);
      nfx[0] >>= 1;

      nf += nfx[1]%2 << 3*i;
      nfx[1] >>= 1;

      nf += nfx[2]%2 << (3*i+2);
      nfx[2] >>= 1;
    }
    nfn[j] = nf;
  }
}
