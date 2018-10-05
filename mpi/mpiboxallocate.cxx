#include "mpi.h"
#include "../include/constants.h"

extern float *xi,*yi,*zi,*xj,*yj,*zj;
extern int *nfn,**nsij,*ncnt,*nek,*irank,*na,*nb;
extern double *tfmm;
extern int *ksdsp,*kscnt,*krdsp,*krcnt;
extern int *nsend,*nrecv;

extern void boxc(int, int, int*);
extern void boxn1(int*, int&, int);
extern void boxn(int, float*, float*, float*);
extern void sort(int);
extern void mpialltoallvi(int*, int*, int*, int*, int*, int*, int);

void boxallocate(int n0, int n1, int n2, int n3) {
  int n,i,nbc,lbi,lbj,lbjd,lbjo,nr,lev;
  int ix,ixp,iy,iyp,iz,izp,jx,jxp,jy,jyp,jz,jzp;
  int ixmin,ixmax,iymin,iymax,izmin,izmax;
  int ic,ie,iv,isrank,jj,idig,ndig,nfdig;
  int nc[3];
  tic = get_time();
  nbmax = 1 << 3*lmax;
  nbned = 1;
  nbnet = 0;

  ksdsp = new int [nprocs];
  kscnt = new int [nprocs];
  krdsp = new int [nprocs];
  krcnt = new int [nprocs];
  nsend = new int [nbmax];
  nrecv = new int [nbmax];
  ncnt = new int [nprocs];
// lbi
  n = n1;
  boxn(n,xi,yi,zi);
  for( i=0; i<n; i++ ) {
    na[i] = nfn[i];
    nb[i] = i+n0;
  }
  sort(n);
  lbi = 0;
  nbc = -1;
  for( i=0; i<n; i++ ) {
    if( na[i] != nbc ) {
      nbc = na[i];
      lbi++;
    }
  }
  nbned = std::max(lbi,nbned);

// lbj 
  n = n3;
  boxn(n,xj,yj,zj);
  for( i=0; i<n; i++ ) {
    na[i] = nfn[i];
    nb[i] = i+n2;
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

  nsij = new int* [nprocs];
  for( i=0; i<nprocs; i++ ) nsij[i] = new int [lbj];
// jcnt
  for( lev=lmax; lev>=2; lev-- ) {
    ixmin = 0;
    ixmax = (1 << lev)-1;
    iymin = 0;
    iymax = (1 << lev)-1;
    izmin = 0;
    izmax = (1 << lev)-1;
    for( i=0; i<nprocs; i++ ) kscnt[i] = 0;
    ndig = 1 << 3*(lmax-lev);
    for( jj=0; jj<lbj; jj++ ) {
      boxc(nfn[jj],3,nc);
      jx = nc[0];
      jy = nc[1];
      jz = nc[2];
      jxp = (jx+2)/2;
      jyp = (jy+2)/2;
      jzp = (jz+2)/2;
      for( i=0; i<nprocs; i++ ) ncnt[i] = 0;
      ncnt[myrank] = 1;
      for( ixp=jxp-1; ixp<=jxp+1; ixp++ ) {
        for( iyp=jyp-1; iyp<=jyp+1; iyp++ ) {
          for( izp=jzp-1; izp<=jzp+1; izp++ ) {
            for( ix=std::max(2*ixp-2,ixmin); ix<=std::min(2*ixp-1,ixmax); ix++ ) {
              for( iy=std::max(2*iyp-2,iymin); iy<=std::min(2*iyp-1,iymax); iy++ ) {
                for( iz=std::max(2*izp-2,izmin); iz<=std::min(2*izp-1,izmax); iz++ ) {
                  nc[0] = ix;
                  nc[1] = iy;
                  nc[2] = iz;
                  boxn1(nc,ie,lev);
                  for( idig=0; idig<ndig; idig++ ) {
                    nfdig = ie*ndig+idig;
                    iv = nek[nfdig];
                    isrank = irank[iv];
                    if( iv != -1 && ncnt[isrank] == 0 ) {
                      nsij[isrank][kscnt[isrank]] = jj;
                      kscnt[isrank]++;
                      ncnt[isrank] = 1;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    ic = 0;
    for( i=0; i<nprocs; i++ ) {
      ksdsp[i] = ic;
      ic += kscnt[i];
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Alltoall(kscnt,1,MPI_INT,krcnt,1,MPI_INT,MPI_COMM_WORLD);
    nr = 0;
    for( i=0; i<nprocs; i++ ) {
      krdsp[i] = nr;
      nr += krcnt[i];
    }
    lbjd = lbj+nr;
    nbned = std::max(lbjd,nbned);
    nbnet += lbjd;
  
  // boxm2m
    lbjo = lbj;
    lbj = 0;
    nbc = -1;
    for( jj=0; jj<lbjo; jj++ ) {
      if( nbc != nfn[jj]/8 ) {
        nbc = nfn[jj]/8;
        nfn[lbj] = nbc;
        lbj++;
      }
    }
  }

  for( i=0; i<nprocs; i++ ) delete[] nsij[i];
  delete[] nsij;

  nbne = nbned;
//  nbne = 1 << 3*lmax;
  nbnet = std::max(nbnet,2*nbne);

  delete[] ksdsp;
  delete[] kscnt;
  delete[] krdsp;
  delete[] krcnt;
  delete[] nsend;
  delete[] nrecv;
  delete[] ncnt;

  toc = get_time();
  tfmm[1] += toc-tic;
}
