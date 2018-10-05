#include "../include/constants.h"

extern float *xj,*yj,*zj;
extern int **ndj,*nej,*nfj,*nfn;
extern double *tfmm;

extern void boxn(int, float*, float*, float*);
extern void sort(int);

void boxdataj(int nj, int& lbj, double& rb) {
  int kl,i,nbc;
  tic = get_time();

  kl = 1 << lmax;
  rb = rd/kl;

  boxn(nj,xj,yj,zj);

  lbj = 0;
  nbc = -1;
  for( i=0; i<nbmax; i++ ) nej[i] = -1;
  for( i=0; i<nj; i++ ) {
    if( nfn[i] != nbc ) {
      nej[nfn[i]] = lbj;
      nfj[lbj] = nfn[i];
      ndj[0][lbj] = i;
      if( lbj > 0 ) ndj[1][lbj-1] = i-1;
      nbc = nfn[i]; 
      lbj++;
    }
  }
  if( lbj > 0 ) ndj[1][lbj-1] = nj-1;

  toc = get_time();
  tfmm[11] += toc-tic;
}
