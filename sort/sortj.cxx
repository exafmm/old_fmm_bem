#include "../include/constants.h"

extern float *xj,*yj,*zj,*gxj,*gyj,*gzj,*vj;
extern int *nbj,*nfn,*na,*nb;
extern double *tfmm;

extern void boxn(int, float*, float*, float*);
extern void sort(int);
extern void sortvar(int, int, float*, int*);
extern void unsortvar(int, int, float*, int*);

void sortj(int& mj) {
  int i;
  tic = get_time();

  boxn(mj,xj,yj,zj);
  for( i=0; i<mj; i++ ) {
    na[i] = nfn[i];
    nb[i] = i;
  }
  sort(mj);
  for( i=0; i<mj; i++ ) {
    nbj[i] = nb[i];
  }

  toc = get_time();
  tfmm[7] += toc-tic;

  sortvar(0,mj,xj,nbj);
  sortvar(0,mj,yj,nbj);
  sortvar(0,mj,zj,nbj);
  sortvar(0,mj,gxj,nbj);
  sortvar(0,mj,gyj,nbj);
  sortvar(0,mj,gzj,nbj);
  sortvar(0,mj,vj,nbj);
}

void unsortj(int& mj) {
  unsortvar(0,mj,xj,nbj);
  unsortvar(0,mj,yj,nbj);
  unsortvar(0,mj,zj,nbj);
  unsortvar(0,mj,gxj,nbj);
  unsortvar(0,mj,gyj,nbj);
  unsortvar(0,mj,gzj,nbj);
  unsortvar(0,mj,vj,nbj);
  tic = get_time();

  toc = get_time();
  tfmm[24] += toc-tic;
}
