#include "../include/constants.h"

extern int *nej,*nfj,*nlbj;
extern double *tfmm;

void boxm2mj(int& lbj, int lev) {
  int lbjo,nbc,j,jj,jb;
  tic = get_time();
  nlbj[lev-1] = nlbj[lev]+lbj;
  lbjo = lbj;
  lbj = 0;
  nbc = -1;
  for( j=0; j<nbmax; j++ ) nej[j] = -1;
  for( jj=0; jj<lbjo; jj++ ) {
    jb = jj+nlbj[lev];
    if( nbc != nfj[jb]/8 ) {
      nbc = nfj[jb]/8;
      nfj[lbj+nlbj[lev-1]] = nbc;
      nej[nbc] = lbj;
      lbj++;
    }
  }

  toc = get_time();
  tfmm[11] += toc-tic;
}
