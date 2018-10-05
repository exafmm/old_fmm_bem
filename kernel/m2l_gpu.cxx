#include "../include/constants.h"

extern int *nfi,*nlbi,*nfj,*nlbj,(*neij)[nebm],*nij,*njb;
extern std::complex<double> (*ax)[mpsym],(*bx)[mpsym];
extern std::complex<double> *ynm,***dnm;
extern float *sr;
extern double *tfmm;
extern int *jbase,*jsize,*istagp,*iendgp,*nvecd,*njj;
extern float *arex,*arey,*arez,*aimx,*aimy,*aimz;
extern float *brex,*brey,*brez,*bimx,*bimy,*bimz;
extern float *ynmre,*ynmim,*dnmre,*dnmim;

extern void boxc(int, int, int*);
extern void boxn1(int*, int&, int);
extern void m2lgpu_(int*, double*, double*, float*, float*, float*, float*,
        float*, float*, float*, float*);

void m2l(int nmp, int mp, int lbi, int lbj, int lev, int ini) {
  int kl,idev,mblok,mpdnm,mmdnm,i,j,m,n,npm,nmm,je,k,nk,nmk,ni,nj,ncall,jj,ii,ib,ij,icall,iblok,jc,jjd;
  int jb,jbd,ix,iy,iz,is,jjdd,jx,jy,jz,isize,im,nc[3];
  double rb,op,tic2;
  tic = get_time();

  tic2 = get_time();
  kl = 1 << lev;
  rb = rd/kl;
  idev = myrank%ngpu;
  mblok = 2*nebm+1;
  mpdnm = (4*mp*mp*mp-mp)/3;
  mmdnm = mpdnm*2*nrbm;

  ynmre = new float [4*mpmax];
  ynmim = new float [4*mpmax];
  dnmre = new float [mmdnm];
  dnmim = new float [mmdnm];

  if( ini != 0 ) {
    for( i=0; i<ini; i++ ) {
      for( j=0; j<nmp; j++ ) {
        ax[i][j] = 0;
      }
    }
  }

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

  ni = 0;
  nj = 0;
  ncall = 0;
  istagp[0] = 0;
  for( jj=0; jj<nbne; jj++ ) njj[jj] = 0;
  for( ii=0; ii<lbi; ii++ ) {
    ni += ((nmp+nblok1)/nblok1)*nblok1;
    if( nij[ii] != 0 ) {
      for( ij=0; ij<nij[ii]; ij++ ) {
        jj = neij[ii][ij];
        if( njj[jj] == 0 ) {
          nj += nmp;
          njj[jj] = 1;
        }
      }
    }
    if ( ni > nimax || nj > njmax ) {
      iendgp[ncall] = ii-1;
      ncall++;
      istagp[ncall] = ii;
      ni = ((nmp+nblok1)/nblok1)*nblok1;
      nj = 0;
      for( jj=0; jj<nbne; jj++ ) njj[jj] = 0;
      for( ij=0; ij<nij[ii]; ij++ ) {
        jj = neij[ii][ij];
        nj += nmp;
        njj[jj] = 1;
      }
    }
  }
  iendgp[ncall] = lbi-1;
  if(lbi != 0) ncall++;
  toc = get_time();
  tfmm[28] += toc-tic2;

  for( icall=0; icall<ncall; icall++ ) {
    tic2 = get_time();
    iblok = 0;
    jc = 0;
    jjd = 0;
    op = 0;
    for( jj=0; jj<nbne; jj++ ) njj[jj] = 0;
    for( ii=istagp[icall]; ii<=iendgp[icall]; ii++ ) {
      for( ij=0; ij<nij[ii]; ij++ ) {
        jj = neij[ii][ij];
        jb = njb[jj];
        if( njj[jj] == 0 ) {
          jbase[jjd] = jc;
          for( j=0; j<nmp; j++ ) {
            brex[jc] = std::real(bx[jb][j]);
            bimx[jc] = std::imag(bx[jb][j]);
            jc++;
          }
          jsize[jjd] = jc-jbase[jjd];
          jjd++;
          njj[jj] = jjd;
        }
      }
      ib = ii+nlbi[lev-1];
      boxc(nfi[ib],3,nc);
      ix = nc[0];
      iy = nc[1];
      iz = nc[2];
      isize = nmp;
      for( is=0; is<isize; is+=nblok1 ) {
        nvecd[iblok*mblok+10] = nij[ii];
        for( ij=0; ij<nij[ii]; ij++ ) {
          jj = neij[ii][ij];
          jbd = jj+nlbj[lev-1];
          jjdd = njj[jj]-1;
          boxc(nfj[jbd],3,nc);
          jx = nc[0];
          jy = nc[1];
          jz = nc[2];
          nc[0] = (ix-jx)+3;
          nc[1] = (iy-jy)+3;
          nc[2] = (iz-jz)+3;
          boxn1(nc,je,3);
          nvecd[iblok*mblok+2*ij+11] = jbase[jjdd];
          nvecd[iblok*mblok+2*ij+12] = je+1;
          op += (double) nblok1*jsize[jjdd];
        }
        iblok++;
      }
    }
    nvecd[0] = idev;
    nvecd[1] = iblok;
    nvecd[2] = mblok;
    nvecd[3] = jc;
    nvecd[4] = 1;
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
      isize = nmp;
      for( is=0; is<isize; is+=nblok1 ) {
        for( i=0; i<std::min(isize-is,nblok1); i++ ) {
          im = iblok*nblok1+i;
          ax[ii][is+i] += std::complex<double>(arex[im],aimx[im]);
        }
        iblok++;
      }
    }
    toc = get_time();
    tfmm[29] += toc-tic2;
  }
  for( jj=0; jj<lbj; jj++ ) {
    jb = njb[jj];
    for( j=0; j<nmp; j++ ) {
      bx[jb][j] = 0;
    }
  }
  delete[] ynmre;
  delete[] ynmim;
  delete[] dnmre;
  delete[] dnmim;

  toc = get_time();
  tfmm[20] += toc-tic;
}
