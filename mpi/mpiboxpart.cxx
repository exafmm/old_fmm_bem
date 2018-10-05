#include "mpi.h"
#include "../include/constants.h"

extern float *xi,*yi,*zi,*xj,*yj,*zj;
extern int *nfn,*nek,*irank,*krdsp,*krcnt,*na,*nb;
extern double *tfmm;

extern void memoryuse();
extern void boxl2g(int&, int&);
extern void boxn(int, float*, float*, float*);
extern void sort(int);

void boxpart(int mi, int mj) {
  int n,i,j,nbc,lbi,lbid,lbj,lbjd,lbk,ii,iwork1,iwork2,ista,iend;
  tic = get_time();

  nbmax = 1 << 3*lmax;
  nek = new int [nbmax];
  irank = new int [nbmax];
  krdsp = new int [nprocs];
  krcnt = new int [nprocs];
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
  lbi = 0;
  nbc = -1;
  for( i=0; i<n; i++ ) {
    if( na[i] != nbc ) {
      nfn[lbi] = na[i];
      nbc = na[i];
      lbi++;
    }
  }

  MPI_Allgather(&lbi,1,MPI_INT,krcnt,1,MPI_INT,MPI_COMM_WORLD);
  lbid = 0;
  for( i=0; i<nprocs; i++ ) {
    krdsp[i] = lbid;
    lbid += krcnt[i];
  }
  MPI_Allgatherv(nfn,lbi,MPI_INT,na,krcnt,krdsp,MPI_INT,MPI_COMM_WORLD);
  sort(lbid);
  nbc = -1;
  for( i=0; i<lbid; i++ ) {
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
  lbj = 0;
  nbc = -1;
  for( i=0; i<n; i++ ) {
    if( na[i] != nbc ) {
      nfn[lbj] = na[i];
      nbc = na[i];
      lbj++;
    }
  }

  MPI_Allgather(&lbj,1,MPI_INT,krcnt,1,MPI_INT,MPI_COMM_WORLD);
  lbjd = 0;
  for( i=0; i<nprocs; i++ ) {
    krdsp[i] = lbjd;
    lbjd += krcnt[i];
  }
  MPI_Allgatherv(nfn,lbj,MPI_INT,na,krcnt,krdsp,MPI_INT,MPI_COMM_WORLD);
  sort(lbjd);
  lbj = 0;
  nbc = -1;
  for( i=0; i<lbjd; i++ ) {
    if( na[i] != nbc ) {
      nbc = na[i];
      if( nek[na[i]] == -1 ) {
        nek[na[i]] = lbk;
        lbk++;
      }
    }
  }

// create partition
  for( i=0; i<nprocs; i++ ) {
    iwork1 = lbk/nprocs;
    iwork2 = lbk%nprocs;
    ista = i*iwork1+std::min(i,iwork2);
    iend = ista+iwork1;
    if( iwork2 > i ) iend++;
    for( j=ista; j<iend; j++ ) {
      irank[j] = i;
    }
  }
 
  delete[] krdsp;
  delete[] krcnt;

  toc = get_time();
  tfmm[4] += toc-tic;
}
