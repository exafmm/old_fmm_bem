#include "../include/constants.h"

extern int *nfi,*nlbi;
extern double *tfmm;

void boxm2mi(int& lbi, int lev) {
  int lbio,nbc,ii,ib;
  tic = get_time();
  nlbi[lev-1] = nlbi[lev]+lbi;
  lbio = lbi;
  lbi = 0;
  nbc = -1;
  for( ii=0; ii<lbio; ii++ ) {
    ib = ii+nlbi[lev];
    if( nbc != nfi[ib]/8 ) {
      nbc = nfi[ib]/8;
      nfi[lbi+nlbi[lev-1]] = nbc;
      lbi++;
    }
  }

  toc = get_time();
  tfmm[10] += toc-tic;
}
