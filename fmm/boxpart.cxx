#include "../include/constants.h"

extern float *xi,*yi,*zi,*xj,*yj,*zj;
extern int *nfn,*nek,*irank,*na,*nb;
extern double *tfmm;

extern void memoryuse();
extern void boxn(int, float*, float*, float*);
extern void sort(int);

void boxpart(int mi, int mj) {
  int n,i,nbc,lbk,ii;
  tic = get_time();

  nbmax = 1 << 3*lmax;
  nek = new int [nbmax];
  irank = new int [nbmax];
  mem = nbmax*2*4;
  memoryuse();

  lbk = 0;
  for( ii=0; ii<nbmax; ii++ ) nek[ii] = -1;

// fill i boxes for nek
  n = mi;
  boxn(n,xi,yi,zi);
  for( i=0; i<n; i++ ) {
    na[i] = nfn[i];
    nb[i] = i;
  }
  sort(n);
  nbc = -1;
  for( i=0; i<n; i++ ) {
    if( na[i] != nbc ) {
      nek[na[i]] = lbk;
      nbc = na[i];
      lbk++;
    }
  }

// fill j boxes for nek
  n = mj;
  boxn(n,xj,yj,zj);
  for( i=0; i<n; i++ ) {
    na[i] = nfn[i];
    nb[i] = i;
  }
  sort(n);
  nbc = -1;
  for( i=0; i<n; i++ ) {
    if( na[i] != nbc ) {
      nbc = na[i];
      if( nek[na[i]] == -1 ) {
        nek[na[i]] = lbk;
        lbk++;
      }
    }
  }

// create partition
  for( i=0; i<lbk; i++ ) {
    irank[i] = 0;
  }

  toc = get_time();
  tfmm[4] += toc-tic;
}
