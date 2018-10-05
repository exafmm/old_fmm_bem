#include "../include/constants.h"

extern int *nfi,*nlbi,*neo;
extern float *anm;
extern std::complex<double> *ynm,***dnm;
extern std::complex<double> (*ax)[mpsym],(*axo)[mpsym];
extern double *tfmm;

extern void boxc(int, int, int*);
extern void boxn1(int*, int&, int);
extern void spharot(std::complex<double>*, std::complex<double>*, int, std::complex<double>**);

void l2l(int nmp, int mp, int lev, int lbi) {
  int kl,lbio,ii,ib,i,nfip,nfic,je,j,k,jk,jks,n,jnk,nk,nks,nc[3];
  double rb,rh;
  std::complex<double> cnm,aaxdd;
  std::complex<double> aax[mpsym],aaxd[mpsym];
  tic = get_time();

  kl = 1 << lev;
  rb = rd/kl;
  lbio = lbi;
  if( lbio < 8 ) lbio = 8;
  for( ii=0; ii<lbio; ii++ ) {
    for( i=0; i<nmp; i++ ) {
      axo[ii][i] = ax[ii][i];
    }
  }
//
  neo = new int [nbmax];
  int nbc;
  nbc = -1;
  lbio = 0;
  for( i=0; i<nbmax; i++ ) neo[i] = -1;
  for( ii=0; ii<lbi; ii++ ) {
    ib = ii+nlbi[lev-1];
    if( nbc != nfi[ib]/8 ) {
      nbc = nfi[ib]/8;
      neo[nbc] = lbio;
      lbio++;
    }
  }
//
  for( ii=0; ii<lbi; ii++ ) {
    ib = ii+nlbi[lev-1];
    nfip = nfi[ib]/8;
    nfic = nfi[ib]%8;
    boxc(nfic,3,nc);
    nc[0] = nc[0]*2+2;
    nc[1] = nc[1]*2+2;
    nc[2] = nc[2]*2+2;
    boxn1(nc,je,3);
    rh = rb*sqrt(3.0)/2;
    ib = neo[nfip];
    for( i=0; i<nmp; i++ ) {
      aaxd[i] = axo[ib][i];
    }
    spharot(aaxd,aax,mp,dnm[je]);
    for( j=0; j<mp; j++ ) {
      for( k=0; k<=j; k++ ) {
        jk = j*j+j+k;
        jks = j*(j+1)/2+k;
        aaxdd = 0;
        for( n=j; n<mp; n++ ) {
          jnk = (n-j)*(n-j)+n-j;
          nk = n*n+n+k;
          nks = n*(n+1)/2+k;
          cnm = anm[jnk]*anm[jk]/anm[nk]*pow(rh,n-j)*ynm[jnk];
          aaxdd += aax[nks]*cnm;
        }
        aaxd[jks] = aaxdd;
      }
    }
    spharot(aaxd,aax,mp,dnm[je+nrbm]);
    for( i=0; i<nmp; i++ ) {
      ax[ii][i] = aax[i];
    }
  }

  delete[] neo;

  toc = get_time();
  tfmm[21] += toc-tic;
}
