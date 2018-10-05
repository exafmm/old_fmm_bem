#include "mpi.h"
#include "../include/constants.h"

extern float *xj,*yj,*zj,*gxj,*gyj,*gzj,*vj;
extern int *nfi,**ndj,*nej,*nfj,*nlbj,**nsij;
extern double *tfmm;
extern int *ksdsp,*kscnt,*krdsp,*krcnt,*lsdsp,*lscnt,*lrdsp,*lrcnt;
extern int *nsend,*nrecv;
extern float *fsend,*frecv;
extern std::complex<float> *csend,*crecv;

extern void memoryuse();
extern void memoryfree();
extern void ijbox(int, int, int, int);
extern void jcnt(int&, int, int, int);
extern void mpialltoallvi(int*, int*, int*, int*, int*, int*, int);
extern void mpisendp2p(int, int, float*);
extern void mpisendm2l(int, int, int, int&, int);

void ifbox(int mi, int mj, int nmp, int lbi, int lbj, int& lbjr, int lev, int neib) {
  int nnmax,nmmax,ii,i,nr,ic,jj,jb,np,lbjd,lbjr2;
  tic = get_time();

  nnmax = nbne*10;
  nmmax = nmp*nnmax;

  ksdsp = new int [nprocs];
  kscnt = new int [nprocs];
  krdsp = new int [nprocs];
  krcnt = new int [nprocs];
  lsdsp = new int [nprocs];
  lscnt = new int [nprocs];
  lrdsp = new int [nprocs];
  lrcnt = new int [nprocs];
  nsend = new int [nnmax];
  nrecv = new int [nnmax];
  csend = new std::complex<float> [nmmax];
  crecv = new std::complex<float> [nmmax];
  mem = nprocs*8*4+nnmax*2*4+npmax*2*4+nmmax*2*8;
  memoryuse();

  toc = get_time();
  tfmm[12] += toc-tic;
  jcnt(nr,lbj,lev,neib);
  tic = get_time();

// nfj
  ic = 0;
  for( ii=0; ii<nprocs; ii++ ) {
    for( i=0; i<kscnt[ii]; i++ ) {
      jj = nsij[ii][i];
      jb = jj+nlbj[lev-1];
      nsend[ic] = nfj[jb];
      ic++;
    }
  }
  mpialltoallvi(nsend,kscnt,ksdsp,nrecv,krcnt,krdsp,nnmax);
  lbjd = lbj;
  for( i=0; i<nr; i++ ) {
    jb = lbjd+nlbj[lev-1];
    nfj[jb] = nrecv[i];
    lbjd++;
  }
// nej
  lbjr = lbj;
  lbjd = lbj;
  for( ii=0; ii<nr; ii++ ) {
    jb = lbjd+nlbj[lev-1];
    if( nej[nfj[jb]] == -1 ) {
      nej[nfj[jb]] = lbjr;
      lbjr++;
    }
    lbjd++;
  }

  if( neib == 0 ) {
// iscnt for xj-vj
    for( i=0; i<nprocs; i++ ) lscnt[i] = 0;
    for( ii=0; ii<nprocs; ii++ ) {
      for( i=0; i<kscnt[ii]; i++ ) {
        jj = nsij[ii][i];
        lscnt[ii] += ndj[1][jj]-ndj[0][jj]+1;
      }
    }
    ic = 0;
    for( i=0; i<nprocs; i++ ) {
      lsdsp[i] = ic;
      ic += lscnt[i];
    }
    MPI_Alltoall(lscnt,1,MPI_INT,lrcnt,1,MPI_INT,MPI_COMM_WORLD);
    np = 0;
    for( i=0; i<nprocs; i++ ) {
      lrdsp[i] = np;
      np += lrcnt[i];
    }
// xj-vj
    toc = get_time();
    tfmm[12] += toc-tic;
    mpisendp2p(mj,np,xj);
    mpisendp2p(mj,np,yj);
    mpisendp2p(mj,np,zj);
    mpisendp2p(mj,np,gxj);
    mpisendp2p(mj,np,gyj);
    mpisendp2p(mj,np,gzj);
    mpisendp2p(mj,np,vj);
    tic = get_time();
// ndj
    ic = 0;
    for( ii=0; ii<nprocs; ii++ ) {
      for( i=0; i<kscnt[ii]; i++ ) {
        jj = nsij[ii][i];
        nsend[ic] = ndj[1][jj]-ndj[0][jj]+1;
        ic++;
      }
    }
    mpialltoallvi(nsend,kscnt,ksdsp,nrecv,krcnt,krdsp,nnmax);
    ic = mj;
    lbjr = lbj;
    for( i=0; i<nr; i++ ) {
      ndj[0][lbjr] = ic;
      ic += nrecv[i];
      ndj[1][lbjr] = ic-1;
      lbjr++;
    }
  } else {
// lscnt for bx-bz (gx-gz)
    for( i=0; i<nprocs; i++ ) {
      lsdsp[i] = nmp*ksdsp[i];
      lscnt[i] = nmp*kscnt[i];
    }
    MPI_Alltoall(lscnt,1,MPI_INT,lrcnt,1,MPI_INT,MPI_COMM_WORLD);
    ic = 0;
    for( i=0; i<nprocs; i++ ) {
      lrdsp[i] = ic;
      ic += lrcnt[i];
    }
    toc = get_time();
    tfmm[12] += toc-tic;
    mpisendm2l(nmp,nr,lbj,lbjr2,lev);
    tic = get_time();
    assert(lbjr==lbjr2);

  }

  toc = get_time();
  tfmm[12] += toc-tic;
  ijbox(lbi,lbjr,lev,neib);
  tic = get_time();

  delete[] ksdsp;
  delete[] kscnt;
  delete[] krdsp;
  delete[] krcnt;
  delete[] lsdsp;
  delete[] lscnt;
  delete[] lrdsp;
  delete[] lrcnt;
  delete[] nsend;
  delete[] nrecv;
  delete[] csend;
  delete[] crecv;
  mem = nprocs*8*4+nnmax*2*4+npmax*2*4+nmmax*2*8;
  memoryfree();

  toc = get_time();
  tfmm[12] += toc-tic;
}
