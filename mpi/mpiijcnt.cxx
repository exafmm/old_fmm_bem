#include "mpi.h"
#include "../include/constants.h"

extern int *nfj,*nlbj,**nsij,*ncnt,*nek,*irank,*ksdsp,*kscnt,*krdsp,*krcnt;
extern double *tfmm;

extern void boxc(int, int, int*);
extern void boxn1(int*, int&, int);

void jcnt(int& nr, int lbj, int lev, int neib) {
  int ix,ixp,iy,iyp,iz,izp,jx,jxp,jy,jyp,jz,jzp;
  int ixmin,ixmax,iymin,iymax,izmin,izmax;
  int i,ic,ie,iv,isrank,jj,jb,idig,ndig,nfdig;
  int nc[3];
  tic = get_time();

  ncnt = new int [nprocs];

  ixmin = 0;
  ixmax = (1 << lev)-1;
  iymin = 0;
  iymax = (1 << lev)-1;
  izmin = 0;
  izmax = (1 << lev)-1;

  for( i=0; i<nprocs; i++ ) kscnt[i] = 0;
  ndig = 1 << 3*(lmax-lev);
  if( neib == 0 ) {
    for( jj=0; jj<lbj; jj++ ) {
      jb = jj+nlbj[lev-1];
      boxc(nfj[jb],3,nc);
      jx = nc[0];
      jy = nc[1];
      jz = nc[2];
      for( i=0; i<nprocs; i++ ) ncnt[i] = 0;
      ncnt[myrank] = 1;
      for( ix=std::max(jx-1,ixmin); ix<=std::min(jx+1,ixmax); ix++ ) {
        for( iy=std::max(jy-1,iymin); iy<=std::min(jy+1,iymax); iy++ ) {
          for( iz=std::max(jz-1,izmin); iz<=std::min(jz+1,izmax); iz++ ) {
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
  } else if( neib == 1 ) {
    for( jj=0; jj<lbj; jj++ ) {
      jb = jj+nlbj[lev-1];
      boxc(nfj[jb],3,nc);
      jx = nc[0];
      jy = nc[1];
      jz = nc[2];
      for( i=0; i<nprocs; i++ ) ncnt[i] = 0;
      ncnt[myrank] = 1;
      for( ie=0; ie<64; ie++ ) {
        boxc(ie,3,nc);
        ix = nc[0];
        iy = nc[1];
        iz = nc[2];
        if( ix < jx-1 || jx+1 < ix || iy < jy-1 || jy+1 < iy || iz < jz-1 || jz+1 < iz ) {
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
  } else if( neib == 2 ) {
    for( jj=0; jj<lbj; jj++ ) {
      jb = jj+nlbj[lev-1];
      boxc(nfj[jb],3,nc);
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
                  if( ix < jx-1 || jx+1 < ix || iy < jy-1 || jy+1 < iy || iz < jz-1 || jz+1 < iz ) {
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

  delete[] ncnt;

  toc = get_time();
  tfmm[14] += toc-tic;
}
