#include "../include/constants.h"

extern int *nej,*nfj,*nlbj;
extern float *anm;
extern std::complex<double> *ynm,***dnm;
extern std::complex<double> (*bx)[mpsym];
extern double *tfmm;

extern void boxc(int, int, int*);
extern void boxn1(int*, int&, int);
extern void spharot(std::complex<double>*, std::complex<double>*, int, std::complex<double>**);

void m2m(int nmp, int mp, int lev, int lbj, int lbjo) {
  int kl,ii,ib,j,jj,nfjp,nfjc,jb,je,k,jk,jks,n,jnk,jnks,nm,nc[3];
  double rb,rh;
  std::complex<double> cnm,bbxdd;
  std::complex<double> bbx[mpsym],bbxd[mpsym];
  tic = get_time();

  kl = 1 << lev;
  rb = rd/kl;
  for( ii=0; ii<lbj; ii++ ) {
    ib = ii+nlbj[lev-1];
    for( j=0; j<nmp; j++ ) {
      bx[ib][j] = 0;
    }
  }
  for( jj=0; jj<lbjo; jj++ ) {
    jb = jj+nlbj[lev];
    nfjp = nfj[jb]/8;
    nfjc = nfj[jb]%8;
    ib = nej[nfjp]+nlbj[lev-1];
    boxc(nfjc,3,nc);
    nc[0] = 4-nc[0]*2;
    nc[1] = 4-nc[1]*2;
    nc[2] = 4-nc[2]*2;
    boxn1(nc,je,3);
    rh = rb*sqrt(3.0)/4;
    for( j=0; j<nmp; j++ ) {
      bbxd[j] = bx[jb][j];
    }
    spharot(bbxd,bbx,mp,dnm[je]);
    for( j=0; j<mp; j++ ) {
      for( k=0; k<=j; k++ ) {
        jk = j*j+j+k;
        jks = j*(j+1)/2+k;
        bbxdd = 0;
        for( n=0; n<=j-abs(k); n++ ) {
          jnk = (j-n)*(j-n)+j-n+k;
          jnks = (j-n)*(j-n+1)/2+k;
          nm = n*n+n;
          cnm = pow(-1.0,n)*anm[nm]*anm[jnk]/anm[jk]*pow(rh,n)*ynm[nm];
          bbxdd += bbx[jnks]*cnm;
        }
        bbxd[jks] = bbxdd;
      }
    }
    spharot(bbxd,bbx,mp,dnm[je+nrbm]);
    for( j=0; j<nmp; j++ ) {
      bx[ib][j] += bbx[j];
    }
  }

  toc = get_time();
  tfmm[19] += toc-tic;
}
