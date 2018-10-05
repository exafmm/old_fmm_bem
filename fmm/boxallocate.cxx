#include "../include/constants.h"

extern float *xi,*yi,*zi,*xj,*yj,*zj;
extern int *nfn,*na,*nb;
extern double *tfmm;

extern void boxn(int, float*, float*, float*);
extern void sort(int);

void boxallocate(int n0, int n1, int n2, int n3) {
  int n,i,nbc,lbi,lbj,lev,nbneit,nbnejt;
  tic = get_time();
  nbmax = 1 << 3*lmax;

  n = n1;
  boxn(n,xi,yi,zi);
  toc = get_time();
  tfmm[29] += toc-tic;
  tic = get_time();
  for( i=0; i<n; i++ ) {
    na[i] = nfn[i];
    nb[i] = i+n0;
  }
  sort(n);
  toc = get_time();
  tfmm[30] += toc-tic;
  tic = get_time();

  lbi = 0;
  nbc = -1;
  for( i=0; i<n; i++ ) {
    if( na[i] != nbc ) {
      lbi++;
      nbc = na[i];
    }
  }
  nbneit = lbi;
  for( lev=lmax-1; lev>=2; lev-- ) {
    nbc = -1;
    for( i=0; i<n; i++ ) {
      if( na[i]/(1 << 3*(lmax-lev)) != nbc ) {
        nbneit++;
        nbc = na[i]/(1 << 3*(lmax-lev));
      }
    }
  }
  toc = get_time();
  tfmm[31] += toc-tic;
  tic = get_time();

  n = n3;
  boxn(n,xj,yj,zj);
  toc = get_time();
  tfmm[29] += toc-tic;
  tic = get_time();
  for( i=0; i<n; i++ ) {
    na[i] = nfn[i];
    nb[i] = i+n2;
  }
  sort(n);
  toc = get_time();
  tfmm[30] += toc-tic;
  tic = get_time();

  lbj = 0;
  nbc = -1;
  for( i=0; i<n; i++ ) {
    if( na[i] != nbc ) {
      lbj++;
      nbc = na[i];
    }
  }
  nbnejt = lbj;
  for( lev=lmax-1; lev>=2; lev-- ) {
    nbc = -1;
    for( i=0; i<n; i++ ) {
      if( na[i]/(1 << 3*(lmax-lev)) != nbc ) {
        nbnejt++;
        nbc = na[i]/(1 << 3*(lmax-lev));
      }
    }
  }
  toc = get_time();
  tfmm[31] += toc-tic;
  tic = get_time();

  nbne = std::max(lbi,lbj);
  nbnet = std::max(nbneit,nbnejt);

  toc = get_time();
  tfmm[1] += toc-tic;
}
