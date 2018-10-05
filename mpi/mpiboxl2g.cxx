#include "mpi.h"
#include "../include/constants.h"

extern int **ndi,*nfi,**ndj,*nej,*nfj,*krdsp,*krcnt,*na,*nb;
extern double *tfmm;

extern void sort(int);

void boxl2g(int& lbi, int& lbj) {
  int i,lbid,nbc,lbjd;
  tic = get_time();

  krdsp = new int [nprocs];
  krcnt = new int [nprocs];

  for( i=0; i<nbmax; i++ ) {
    nej[i] = 0;
  }

  MPI_Allgather(&lbi,1,MPI_INT,krcnt,1,MPI_INT,MPI_COMM_WORLD);
  lbid = 0;
  for( i=0; i<nprocs; i++ ) {
    krdsp[i] = lbid;
    lbid += krcnt[i];
  }

  MPI_Allgatherv(nfi,lbi,MPI_INT,na,krcnt,krdsp,MPI_INT,MPI_COMM_WORLD);
  sort(lbid);

  lbi = 0;
  nbc = -1;
  for( i=0; i<lbid; i++ ) {
    if( na[i] != nbc ) {
      nfi[lbi] = na[i];
      nbc = na[i];
      lbi++;
    }
  }

  MPI_Allgather(&lbj,1,MPI_INT,krcnt,1,MPI_INT,MPI_COMM_WORLD);
  lbjd = 0;
  for( i=0; i<nprocs; i++ ) {
    krdsp[i] = lbjd;
    lbjd += krcnt[i];
  }
  MPI_Allgatherv(nfj,lbj,MPI_INT,na,krcnt,krdsp,MPI_INT,MPI_COMM_WORLD);
  sort(lbjd);
  lbj = 0;
  nbc = -1;
  for( i=0; i<lbjd; i++ ) {
    if( na[i] != nbc ) {
      nej[na[i]] = lbj;
      nfj[lbj] = na[i];
      nbc = na[i];
      lbj++;
    }
  }

  delete[] krdsp;
  delete[] krcnt;

  toc = get_time();
  tfmm[5] += toc-tic;
}
