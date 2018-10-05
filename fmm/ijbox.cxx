#include "../include/constants.h"

extern int *nei,*nfi,*nlbi,*nej,*nfj,*nlbj,(*neij)[nebm],*nij,*njb;
extern double *tfmm;

extern void boxc(int, int, int*);
extern void boxn1(int*, int&, int);

void ijbox(int lbi, int lbj, int lev, int neib) {
  int jxmin,jxmax,jymin,jymax,jzmin,jzmax,ii,ib,jj,jb,ix,iy,iz,jx,jy,jz,je;
  int ixp,iyp,izp,jxp,jyp,jzp,nc[3];
  tic = get_time();

  jxmin = 1000000;
  jxmax = -1000000;
  jymin = 1000000;
  jymax = -1000000;
  jzmin = 1000000;
  jzmax = -1000000;
  for( jj=0; jj<lbj; jj++ ) {
    jb = jj+nlbj[lev-1];
    boxc(nfj[jb],3,nc);
    jxmin = std::min(jxmin,nc[0]);
    jxmax = std::max(jxmax,nc[0]);
    jymin = std::min(jymin,nc[1]);
    jymax = std::max(jymax,nc[1]);
    jzmin = std::min(jzmin,nc[2]);
    jzmax = std::max(jzmax,nc[2]);
  }
  if( neib == 0 ) {
    for( ii=0; ii<lbi; ii++ ) {
      ib = ii+nlbi[lev-1];
      nij[ii] = 0;
      boxc(nfi[ib],3,nc);
      ix = nc[0];
      iy = nc[1];
      iz = nc[2];
      for( jx=std::max(ix-1,jxmin); jx<=std::min(ix+1,jxmax); jx++ ) {
        for( jy=std::max(iy-1,jymin); jy<=std::min(iy+1,jymax); jy++ ) {
          for( jz=std::max(iz-1,jzmin); jz<=std::min(iz+1,jzmax); jz++ ) {
            nc[0] = jx;
            nc[1] = jy;
            nc[2] = jz;
            boxn1(nc,je,lev);
            jj = nej[je];
            if( jj != -1 ) {
              neij[ii][nij[ii]] = jj;
              nij[ii]++;
            }
          }
        }
      }
    }
  } else if( neib == 1 ) {
    for( ii=0; ii<lbi; ii++ ) {
      ib = ii+nlbi[lev-1];
      nij[ii] = 0;
      boxc(nfi[ib],3,nc);
      ix = nc[0];
      iy = nc[1];
      iz = nc[2];
      for( jj=0; jj<lbj; jj++ ) {
        jb = jj+nlbj[lev-1];
        boxc(nfj[jb],3,nc);
        jx = nc[0];
        jy = nc[1];
        jz = nc[2];
        if( jx < ix-1 || ix+1 < jx || jy < iy-1 || iy+1 < jy || jz < iz-1 || iz+1 < jz ) {
          neij[ii][nij[ii]] = jj;
          nij[ii]++;
        }
      }
    }
  } else if( neib == 2 ) {
    for( ii=0; ii<lbi; ii++ ) {
      ib = ii+nlbi[lev-1];
      nij[ii] = 0;
      boxc(nfi[ib],3,nc);
      ix = nc[0];
      iy = nc[1];
      iz = nc[2];
      ixp = (ix+2)/2;
      iyp = (iy+2)/2;
      izp = (iz+2)/2;
      for( jxp=ixp-1; jxp<=ixp+1; jxp++ ) {
        for( jyp=iyp-1; jyp<=iyp+1; jyp++ ) {
          for( jzp=izp-1; jzp<=izp+1; jzp++ ) {
            for( jx=std::max(2*jxp-2,jxmin); jx<=std::min(2*jxp-1,jxmax); jx++ ) {
              for( jy=std::max(2*jyp-2,jymin); jy<=std::min(2*jyp-1,jymax); jy++ ) {
                for( jz=std::max(2*jzp-2,jzmin); jz<=std::min(2*jzp-1,jzmax); jz++ ) {
                  if( jx < ix-1 || ix+1 < jx || jy < iy-1 || iy+1 < jy || jz < iz-1 || iz+1 < jz ) {
                    nc[0] = jx;
                    nc[1] = jy;
                    nc[2] = jz;
                    boxn1(nc,je,lev);
                    jj = nej[je];
                    if( jj != -1 ) {
                      neij[ii][nij[ii]] = jj;
                      nij[ii]++;
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
  for( jj=0; jj<lbj; jj++ ) {
    jb = jj+nlbj[lev-1];
    njb[jj] = nej[nfj[jb]]+nlbj[lev-1];
  }

  toc = get_time();
  tfmm[13] += toc-tic;
}
