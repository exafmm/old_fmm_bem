#include "../include/constants.h"

extern int *nej,*nfj,*nlbj;
extern std::complex<double> (*bx)[mpsym];
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

void m2m(int nmp, int mp, int lev, int lbj, int lbjo) {
  int kl,idev,mblok,mpdnm,mmdnm,m,n,npm,nmm,je,k,nk,nmk,ii,ib,j,nj,ncall,jj,icall;
  int iblok,jc,nfjc,jb,jbase,jsize,nfjp,jm,nc[3];
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

  for( ii=0; ii<lbj; ii++ ) {
    ib = ii+nlbj[lev-1];
    for( j=0; j<nmp; j++ ) {
      bx[ib][j] = 0;
    }
  }

  nj = 0;
  ncall = 0;
  istagp[0] = 0;
  for( jj=0; jj<lbjo; jj++ ) {
    nj += nblok1;
    if( nj > njmax ) {
      iendgp[ncall] = jj-1;
      ncall++;
      istagp[ncall] = jj;
      nj = nblok1;
    }
  }
  iendgp[ncall] = lbjo-1;
  if(lbjo != 0) ncall++;
  toc = get_time();
  tfmm[28] += toc-tic2;

  for( icall=0; icall<ncall; icall++ ) {
    tic2 = get_time();
    iblok = 0;
    jc = 0;
    op = 0;
    for( jj=istagp[icall]; jj<=iendgp[icall]; jj++ ) {
      jb = jj+nlbj[lev];
      nfjc = nfj[jb]%8;
      boxc(nfjc,3,nc);
      nc[0] = 4-nc[0]*2;
      nc[1] = 4-nc[1]*2;
      nc[2] = 4-nc[2]*2;
      boxn1(nc,je,3);
      jbase = jc;
      for( j=0; j<nmp; j++ ) {
        brex[jc] = std::real(bx[jb][j]);
        bimx[jc] = std::imag(bx[jb][j]);
        jc++;
      }
      jsize = jc-jbase;
      nvecd[iblok*mblok+10] = 1;
      nvecd[iblok*mblok+11] = jbase;
      nvecd[iblok*mblok+12] = je+1;
      op += nblok1*jsize;
      iblok++;
    }

    nvecd[0] = idev;
    nvecd[1] = iblok;
    nvecd[2] = mblok;
    nvecd[3] = jc;
    nvecd[4] = 0;
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
    for( jj=istagp[icall]; jj<=iendgp[icall]; jj++ ) {
      jb = jj+nlbj[lev];
      nfjp = nfj[jb]/8;
      ib = nej[nfjp]+nlbj[lev-1];
      for( j=0; j<nmp; j++ ) {
        jm = iblok*nblok1+j;
        bx[ib][j] += std::complex<double>(arex[jm],aimx[jm]);
      }
      iblok++;
    }
    toc = get_time();
    tfmm[29] += toc-tic2;
  }
  delete[] ynmre;
  delete[] ynmim;
  delete[] dnmre;
  delete[] dnmim;

  toc = get_time();
  tfmm[19] += toc-tic;
}
