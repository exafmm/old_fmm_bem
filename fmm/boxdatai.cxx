#include "../include/constants.h"

extern float *xi,*yi,*zi;
extern int **ndi,*nfi,*nfn;
extern double *tfmm;

extern void boxn(int, float*, float*, float*);
extern void sort(int);

void boxdatai(int ni, int& lbi, double& rb) {
  int kl,i,nbc;
  tic = get_time();

  kl = 1 << lmax;
  rb = rd/kl;

  boxn(ni,xi,yi,zi);

  lbi = 0;
  nbc = -1;
  for( i=0; i<ni; i++ ) {
    if( nfn[i] != nbc ) {
      nfi[lbi] = nfn[i];
      ndi[0][lbi] = i;
      if( lbi > 0 ) ndi[1][lbi-1] = i-1;
      nbc = nfn[i]; 
      lbi++;
    }
  }
  if( lbi > 0 ) ndi[1][lbi-1] = ni-1;

  toc = get_time();
  tfmm[10] += toc-tic;
}
