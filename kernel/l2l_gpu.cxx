#include "../include/constants.h"

extern int *nfi,*nlbi,*neo;
extern std::complex<double> (*ax)[mpsym],(*axo)[mpsym];
extern std::complex<double> *ynm,***dnm;
extern double *tfmm;
extern int *istagp,*iendgp,*nvecd;
extern float *arex,*aimx;
extern float *brex,*bimx;
extern float *ynmre,*ynmim,*dnmre,*dnmim;

extern void boxc(int, int, int*);
extern void boxn1(int*, int&, int);
extern void m2lgpu_(int*, double*, double*, float*, float*, float*, float*,
        float*, float*, float*, float*);

void l2l(int nmp, int mp, int lev, int lbi) {
  int kl,idev,mblok,mpdnm,mmdnm,m,n,npm,nmm,je,k,nk,nmk,lbio,ii,i,ni;
  int ncall,icall,iblok,ic,nfip,nfic,ib,jbase,jsize,im,nc[3];
  double rb,op,tic2;
  tic = get_time();

  tic2 = get_time();
  kl = 1 << lev;
  rb = rd/kl;
  idev = myrank%ngpu;
  mblok = 3;
  mpdnm = (4*mp*mp*mp-mp)/3;
  mmdnm = mpdnm*2*nrbm;

  ynmre = new float [4*mpmax];
  ynmim = new float [4*mpmax];
  dnmre = new float [mmdnm];
  dnmim = new float [mmdnm];

  for( m=0; m<2*mp; m++ ) {
    for( n=m; n<2*mp; n++ ) {
      npm = n*n+n+m;
      nmm = n*n+n-m;
      ynmre[npm] = std::real(ynm[npm]);
      ynmre[nmm] = std::real(ynm[nmm]);
      ynmim[npm] = std::imag(ynm[npm]);
      ynmim[nmm] = std::imag(ynm[nmm]);
    }
  }
  for( je=0; je<2*nrbm; je++ ) {
    for( n=0; n<mp; n++ ) {
      for( m=0; m<=n; m++ ) {
        for( k=-n; k<=n; k++ ) {
          nk = n*(n+1)+k;
          nmk = (4*n*n*n+6*n*n+5*n)/3+m*(2*n+1)+k+je*mpdnm;
          dnmre[nmk] = std::real(dnm[je][m][nk]);
          dnmim[nmk] = std::imag(dnm[je][m][nk]);
        }
      }
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
  lbio = lbi;
  if( lbio < 8 ) lbio = 8;
  for( ii=0; ii<lbio; ii++ ) {
    for( i=0; i<nmp; i++ ) {
      axo[ii][i] = ax[ii][i];
    }
  }

  ni = 0;
  ncall = 0;
  istagp[0] = 0;
  for( ii=0; ii<lbi; ii++ ) {
    ni += nblok1;
    if( ni > njmax ) {
      iendgp[ncall] = ii-1;
      ncall++;
      istagp[ncall] = ii;
      ni = nblok1;
    }
  }
  iendgp[ncall] = lbi-1;
  if(lbi != 0) ncall++;
  toc = get_time();
  tfmm[28] += toc-tic2;

  for( icall=0; icall<ncall; icall++ ) {
    iblok = 0;
    ic = 0;
    op = 0;
    for( ii=istagp[icall]; ii<=iendgp[icall]; ii++ ) {
      ib = ii+nlbi[lev-1];
      nfip = nfi[ib]/8;
      nfic = nfi[ib]%8;
      boxc(nfic,3,nc);
      nc[0] = nc[0]*2+2;
      nc[1] = nc[1]*2+2;
      nc[2] = nc[2]*2+2;
      boxn1(nc,je,3);
      ib = neo[nfip];
      jbase = ic;
      for( i=0; i<nmp; i++ ) {
        brex[ic] = std::real(axo[ib][i]);
        bimx[ic] = std::imag(axo[ib][i]);
        ic++;
      }
      jsize = ic-jbase;
      nvecd[iblok*mblok+10] = 1;
      nvecd[iblok*mblok+11] = jbase;
      nvecd[iblok*mblok+12] = je+1;
      op += nblok1*jsize;
      iblok++;
    }

    nvecd[0] = idev;
    nvecd[1] = iblok;
    nvecd[2] = mblok;
    nvecd[3] = ic;
    nvecd[4] = 2;
    nvecd[5] = myrank;
    nvecd[6] = mp;
    nvecd[7] = nrbm;
    toc = get_time();
    tfmm[29] += toc-tic2;
    tic2 = get_time();
    m2lgpu_(nvecd,&op,&rb,arex,aimx,brex,bimx,ynmre,ynmim,dnmre,dnmim);
    toc = get_time();
    tfmm[30] += toc-tic2;
    tic2 = get_time();
    iblok = 0;
    for( ii=istagp[icall]; ii<=iendgp[icall]; ii++ ) {
      for( i=0; i<nmp; i++ ) {
        im = iblok*nblok1+i;
        ax[ii][i] = std::complex<double>(arex[im],aimx[im]);
      }
      iblok++;
    }
    toc = get_time();
    tfmm[29] += toc-tic2;
  }
  delete[] neo;
  delete[] ynmre;
  delete[] ynmim;
  delete[] dnmre;
  delete[] dnmim;

  toc = get_time();
  tfmm[21] += toc-tic;
}
