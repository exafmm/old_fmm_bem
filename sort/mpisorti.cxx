#include "mpi.h"
#include "../include/constants.h"

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi,*vi,*vo;
extern int *nbi,*nfn,*nek,*irank,*isort,*isdsp,*iscnt,*irdsp,*ircnt,*na,*nb;
extern double *tfmm;

extern void boxn(int, float*, float*, float*);
extern void sort(int);
extern void sortvar(int, int, float*, int*);
extern void unsortvar(int, int, float*, int*);
extern void mpisortvar(int*, int*, int*, int*, float*);
extern void mpiunsortvar(int*, int*, int*, int*, float*);

void sorti(int& mi) {
  int i,iv,lbi,nbc,ista=0,ic;
  tic = get_time();

  boxn(mi,xi,yi,zi);
  for( i=0; i<mi; i++ ) {
    iv = nek[nfn[i]];
    na[i] = irank[iv];
    nb[i] = i;
  }
  sort(mi);
  for( i=0; i<mi; i++ ) {
    isort[i] = nb[i];
  }

  lbi = 0;
  nbc = -1;
  for( i=0; i<nprocs; i++ ) iscnt[i] = 0;
  for( i=0; i<mi; i++ ) {
    if( na[i] != nbc ) {
      lbi++;
      if( lbi >= 2 ) {
        iscnt[na[i-1]] = i-ista;
      }
      ista = i;
      nbc = na[i];
    }
  }
  iscnt[na[mi-1]] = mi-ista;
  ic = 0;
  for( i=0; i<nprocs; i++ ) {
    isdsp[i] = ic;
    ic += iscnt[i];
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Alltoall(iscnt,1,MPI_INT,ircnt,1,MPI_INT,MPI_COMM_WORLD);

  toc = get_time();
  tfmm[6] += toc-tic;
  sortvar(0,mi,xi,isort);
  sortvar(0,mi,yi,isort);
  sortvar(0,mi,zi,isort);
  sortvar(0,mi,gxi,isort);
  sortvar(0,mi,gyi,isort);
  sortvar(0,mi,gzi,isort);
  tic = get_time();

  mi = 0;
  for( i=0; i<nprocs; i++ ) {
    irdsp[i] = mi;
    mi += ircnt[i];
  }

  toc = get_time();
  tfmm[6] += toc-tic;
  mpisortvar(isdsp,iscnt,irdsp,ircnt,xi);
  mpisortvar(isdsp,iscnt,irdsp,ircnt,yi);
  mpisortvar(isdsp,iscnt,irdsp,ircnt,zi);
  mpisortvar(isdsp,iscnt,irdsp,ircnt,gxi);
  mpisortvar(isdsp,iscnt,irdsp,ircnt,gyi);
  mpisortvar(isdsp,iscnt,irdsp,ircnt,gzi);
  tic = get_time();

  boxn(mi,xi,yi,zi);
  for( i=0; i<mi; i++ ) {
    na[i] = nfn[i];
    nb[i] = i;
  }
  sort(mi);
  for( i=0; i<mi; i++ ) nbi[i] = nb[i];

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
  int i;

  unsortvar(0,mi,xi,nbi);
  unsortvar(0,mi,yi,nbi);
  unsortvar(0,mi,zi,nbi);
  unsortvar(0,mi,gxi,nbi);
  unsortvar(0,mi,gyi,nbi);
  unsortvar(0,mi,gzi,nbi);
  unsortvar(0,mi,vi,nbi);

  mpiunsortvar(isdsp,iscnt,irdsp,ircnt,xi);
  mpiunsortvar(isdsp,iscnt,irdsp,ircnt,yi);
  mpiunsortvar(isdsp,iscnt,irdsp,ircnt,zi);
  mpiunsortvar(isdsp,iscnt,irdsp,ircnt,gxi);
  mpiunsortvar(isdsp,iscnt,irdsp,ircnt,gyi);
  mpiunsortvar(isdsp,iscnt,irdsp,ircnt,gzi);
  mpiunsortvar(isdsp,iscnt,irdsp,ircnt,vi);
  tic = get_time();

  mi = 0;
  for( i=0; i<nprocs; i++ ) mi += iscnt[i];

  toc = get_time();
  tfmm[23] += toc-tic;
  unsortvar(0,mi,xi,isort);
  unsortvar(0,mi,yi,isort);
  unsortvar(0,mi,zi,isort);
  unsortvar(0,mi,gxi,isort);
  unsortvar(0,mi,gyi,isort);
  unsortvar(0,mi,gzi,isort);
  unsortvar(0,mi,vi,isort);
}
