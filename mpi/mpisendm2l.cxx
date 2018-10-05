#include "mpi.h"
#include "../include/constants.h"

extern int *nej,*nfj,*nlbj,*nfo,*nij,**nsij;
extern std::complex<double> (*bx)[mpsym];
extern double *tfmm;
extern int *kscnt,*lsdsp,*lscnt,*lrdsp,*lrcnt;
extern std::complex<float> *csend,*crecv;

extern void mpialltoallvc(std::complex<float>*, int*, int*, std::complex<float>*, int*, int*, int);

void mpisendm2l(int nmp, int nr, int lbj, int& lbjr, int lev) {
  int nnmax,nmmax,ic,jj,j,jb,jbd,k,ii,i;
  tic = get_time();

  nnmax = nbne*nprocs;
  nmmax = nmp*nnmax;

// bx
  ic = 0;
  for( ii=0; ii<nprocs; ii++ ) {
    for( i=0; i<kscnt[ii]; i++ ) {
      jbd = nsij[ii][i]+nlbj[lev-1];
      jb = nej[nfj[jbd]]+nlbj[lev-1];
      for( k=0; k<nmp; k++ ) {
        csend[ic] = bx[jb][k];
        ic++;
      }
    }
  }
  mpialltoallvc(csend,lscnt,lsdsp,crecv,lrcnt,lrdsp,nmmax);
  for( i=0; i<nbne; i++ ) nij[i] = 0;
  for( jj=0; jj<lbj; jj++ ) {
    nij[jj]++;
  }
  jj = lbj;
  lbjr = lbj;
  ic = 0;
  for( i=0; i<nr; i++ ) {
    jbd = jj+nlbj[lev-1];
    j = nej[nfj[jbd]];
    jb = nej[nfj[jbd]]+nlbj[lev-1];
    for( k=0; k<nmp; k++ ) {
      bx[jb][k] += crecv[ic];
      ic++;
    }
    if( nij[j] == 0 ) {
      nfo[lbjr] = nfj[jbd];
      nij[j]++;
      lbjr++;
    }
    jj++;
  }
  for( jj=lbj; jj<lbjr; jj++ ) {
    jbd = jj+nlbj[lev-1];
    nfj[jbd] = nfo[jj];
  }

  toc = get_time();
  tfmm[16] += toc-tic;
}
