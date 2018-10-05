#include "../include/constants.h"

extern int *nfi,*nlbi,*nfj,*nlbj,(*neij)[nebm],*nij,*njb;
extern float *sr;
extern std::complex<double> (*ax)[mpsym],(*bx)[mpsym];
extern std::complex<double> *ynm,***dnm;
extern double *tfmm;

extern void boxc(int, int, int*);
extern void boxn1(int*, int&, int);
extern void spharot(std::complex<double>*, std::complex<double>*, int, std::complex<double>**);

void m2l(int nmp, int mp, int lbi, int lbj, int lev, int ini) {
  int kl,i,j,ii,ib,ix,iy,iz,ij,jj,jb,jbd,jx,jy,jz,je,k,jk,jks,n,nk,nks,jkn,jnk,nc[3];
  double rb,xijc,yijc,zijc,rh,rhj,rhjk,rhjn;
  std::complex<double> aax[mpsym],bbx[mpsym];
  std::complex<double> aaxd[mpsym],bbxd[mpsym];
  std::complex<double> cnm,aaxdd;
  tic = get_time();

  kl = 1 << lev;
  rb = rd/kl;
  if( ini != 0 ) {
    for( i=0; i<ini; i++ ) {
      for( j=0; j<nmp; j++ ) {
        ax[i][j] = 0;
      }
    }
  }
  for( ii=0; ii<lbi; ii++ ) {
    ib = ii+nlbi[lev-1];
    boxc(nfi[ib],3,nc);
    ix = nc[0];
    iy = nc[1];
    iz = nc[2];
    for( ij=0; ij<nij[ii]; ij++ ) {
      jj = neij[ii][ij];
      jb = njb[jj];
      jbd = jj+nlbj[lev-1];
      for( j=0; j<nmp; j++ ) {
        bbx[j] = bx[jb][j];
      }
      boxc(nfj[jbd],3,nc);
      jx = nc[0];
      jy = nc[1];
      jz = nc[2];
      xijc = (ix-jx)*rb;
      yijc = (iy-jy)*rb;
      zijc = (iz-jz)*rb;
      nc[0] = (ix-jx)+3;
      nc[1] = (iy-jy)+3;
      nc[2] = (iz-jz)+3;
      boxn1(nc,je,3);
      rh = sqrt(xijc*xijc+yijc*yijc+zijc*zijc)+eps;
      spharot(bbx,bbxd,mp,dnm[je]);
      rhj = 1;
      for( j=0; j<mp; j++ ) {
        rhjk = rhj;
        rhj *= rh;
        for( k=0; k<=j; k++ ) {
          jk = j*j+j+k;
          jks = j*(j+1)/2+k;
          aaxdd = 0;
          rhjn = rhjk;
          rhjk *= rh;
          for( n=abs(k); n<mp; n++ ) {
            rhjn *= rh;
            nk = n*n+n+k;
            nks = n*(n+1)/2+k;
            jkn = jk*mp*mp+nk;
            jnk = (j+n)*(j+n)+j+n;
            cnm = sr[jkn]/rhjn*ynm[jnk];
            aaxdd += bbxd[nks]*cnm;
          }
          aaxd[jks] = aaxdd;
        }
      }
      spharot(aaxd,aax,mp,dnm[je+nrbm]);
      for( j=0; j<nmp; j++ ) {
        ax[ii][j] += aax[j];
      }
    }
  }
  for( jj=0; jj<lbj; jj++ ) {
    jb = njb[jj];
    for( j=0; j<nmp; j++ ) {
      bx[jb][j] = 0;
    }
  }

  toc = get_time();
  tfmm[20] += toc-tic;
}
