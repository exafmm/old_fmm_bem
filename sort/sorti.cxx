#include "../include/constants.h"

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi,*vi;
extern int *nbi,*nfn,*na,*nb;
extern double *tfmm;

extern void boxn(int, float*, float*, float*);
extern void sort(int);
extern void sortvar(int, int, float*, int*);
extern void unsortvar(int, int, float*, int*);

void sorti(int& mi) {
  int i;
  tic = get_time();

  boxn(mi,xi,yi,zi);
  for( i=0; i<mi; i++ ) {
    na[i] = nfn[i];
    nb[i] = i;
  }
  sort(mi);
  for( i=0; i<mi; i++ ) {
    nbi[i] = nb[i];
  }

  toc = get_time();
  tfmm[6] += toc-tic;

  sortvar(0,mi,xi,nbi);
  sortvar(0,mi,yi,nbi);
  sortvar(0,mi,zi,nbi);
  sortvar(0,mi,gxi,nbi);
  sortvar(0,mi,gyi,nbi);
  sortvar(0,mi,gzi,nbi);
}

void unsorti(int& mi) {
  unsortvar(0,mi,xi,nbi);
  unsortvar(0,mi,yi,nbi);
  unsortvar(0,mi,zi,nbi);
  unsortvar(0,mi,gxi,nbi);
  unsortvar(0,mi,gyi,nbi);
  unsortvar(0,mi,gzi,nbi);
  unsortvar(0,mi,vi,nbi);
  tic = get_time();

  toc = get_time();
  tfmm[23] += toc-tic;
}
