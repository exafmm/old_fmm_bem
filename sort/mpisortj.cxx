#include "mpi.h"
#include "../include/constants.h"

extern float *xj,*yj,*zj,*gxj,*gyj,*gzj,*vj;
extern int *nbj,*nfn,*nek,*irank,*jsort,*jsdsp,*jscnt,*jrdsp,*jrcnt,*na,*nb;
extern double *tfmm;

extern void boxn(int, float*, float*, float*);
extern void sort(int);
extern void sortvar(int, int, float*, int*);
extern void unsortvar(int, int, float*, int*);
extern void mpisortvar(int*, int*, int*, int*, float*);
extern void mpiunsortvar(int*, int*, int*, int*, float*);

void sortj(int& mj) {
  int i,jv,lbj,nbc,ista=0,ic;
  tic = get_time();

  boxn(mj,xj,yj,zj);
  for( i=0; i<mj; i++ ) {
    jv = nek[nfn[i]];
    na[i] = irank[jv];
    nb[i] = i;
  }
  sort(mj);
  for( i=0; i<mj; i++ ) {
    jsort[i] = nb[i];
  }

  lbj = 0;
  nbc = -1;
  for( i=0; i<nprocs; i++ ) jscnt[i] = 0;
  for( i=0; i<mj; i++ ) {
    if( na[i] != nbc ) {
      lbj++;
      if( lbj >= 2 ) {
        jscnt[na[i-1]] = i-ista;
      }
      ista = i;
      nbc = na[i];
    }
  }
  jscnt[na[mj-1]] = mj-ista;
  ic = 0;
  for( i=0; i<nprocs; i++ ) {
    jsdsp[i] = ic;
    ic += jscnt[i];
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Alltoall(jscnt,1,MPI_INT,jrcnt,1,MPI_INT,MPI_COMM_WORLD);

  toc = get_time();
  tfmm[7] += toc-tic;
  sortvar(0,mj,xj,jsort);
  sortvar(0,mj,yj,jsort);
  sortvar(0,mj,zj,jsort);
  sortvar(0,mj,gxj,jsort);
  sortvar(0,mj,gyj,jsort);
  sortvar(0,mj,gzj,jsort);
  sortvar(0,mj,vj,jsort);
  tic = get_time();

  mj = 0;
  for( i=0; i<nprocs; i++ ) {
    jrdsp[i] = mj;
    mj += jrcnt[i];
  }
  if(mj==0) std::cout << myrank << std::endl;

  toc = get_time();
  tfmm[7] += toc-tic;
  mpisortvar(jsdsp,jscnt,jrdsp,jrcnt,xj);
  mpisortvar(jsdsp,jscnt,jrdsp,jrcnt,yj);
  mpisortvar(jsdsp,jscnt,jrdsp,jrcnt,zj);
  mpisortvar(jsdsp,jscnt,jrdsp,jrcnt,gxj);
  mpisortvar(jsdsp,jscnt,jrdsp,jrcnt,gyj);
  mpisortvar(jsdsp,jscnt,jrdsp,jrcnt,gzj);
  mpisortvar(jsdsp,jscnt,jrdsp,jrcnt,vj);
  tic = get_time();

  boxn(mj,xj,yj,zj);
  for( i=0; i<mj; i++ ) {
    na[i] = nfn[i];
    nb[i] = i;
  }
  sort(mj);
  for( i=0; i<mj; i++ ) nbj[i] = nb[i];

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
  int i;

  unsortvar(0,mj,xj,nbj);
  unsortvar(0,mj,yj,nbj);
  unsortvar(0,mj,zj,nbj);
  unsortvar(0,mj,gxj,nbj);
  unsortvar(0,mj,gyj,nbj);
  unsortvar(0,mj,gzj,nbj);
  unsortvar(0,mj,vj,nbj);

  mpiunsortvar(jsdsp,jscnt,jrdsp,jrcnt,xj);
  mpiunsortvar(jsdsp,jscnt,jrdsp,jrcnt,yj);
  mpiunsortvar(jsdsp,jscnt,jrdsp,jrcnt,zj);
  mpiunsortvar(jsdsp,jscnt,jrdsp,jrcnt,gxj);
  mpiunsortvar(jsdsp,jscnt,jrdsp,jrcnt,gyj);
  mpiunsortvar(jsdsp,jscnt,jrdsp,jrcnt,gzj);
  mpiunsortvar(jsdsp,jscnt,jrdsp,jrcnt,vj);
  tic = get_time();

  mj = 0;
  for( i=0; i<nprocs; i++ ) mj += jscnt[i];

  toc = get_time();
  tfmm[24] += toc-tic;
  unsortvar(0,mj,xj,jsort);
  unsortvar(0,mj,yj,jsort);
  unsortvar(0,mj,zj,jsort);
  unsortvar(0,mj,gxj,jsort);
  unsortvar(0,mj,gyj,jsort);
  unsortvar(0,mj,gzj,jsort);
  unsortvar(0,mj,vj,jsort);
}
