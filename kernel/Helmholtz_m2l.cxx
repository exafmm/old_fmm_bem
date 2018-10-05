#include "../include/constants.h"

extern int *nfi,*nlbi,*nej,*nfj,*nlbj,(*neij)[nebm],*nij,*njb;
extern float *sr,(*ytop)[mp+1],(*m2e)[mp+1][mp+1];
extern double (*rdpi2)[mp+1][2*mp+1];
extern double (*rdmpi2)[mp+1][2*mp+1];
extern double (*rdsq3)[mp+1][2*mp+1];
extern double (*rdmsq3)[mp+1][2*mp+1];
extern std::complex<double> (*ax)[mpsym],(*bx)[mpsym];
extern std::complex<double> *ynm,***dnm;
extern int *iexp1,*iexp2;
extern std::complex<double> *fexpe,*fexpo,*fexpb;
extern std::complex<double> (*xs)[3],(*ys)[3];
extern std::complex<double> **lexp1,**lexp2;
extern double (*zs)[3];
extern double *tfmm;

extern void boxc(int, int, int*);
extern void boxn1(int*, int&, int);
extern void spharot(std::complex<double>*, std::complex<double>*, int, std::complex<double>**);
extern void precalch(int,double*,double*,int*,int*);

void Helmholtz_m2l(int nmp, int mp, int lbi, int lbj, int lev, int ini) {
  int kl,i,j,ii,ib,jj,jb,jx,jy,jz,je,l,m,n,nm,ntot,ic,nbox;
  int ival,nftot,nptot,nexte,nexto,next,mm,islast;
  int jxmin,jxmax,jymin,jymax,jzmin,jzmax;
  int ix,iy,iz,ixp,iyp,izp,jxp,jyp,jzp,nc[3];
  int nfour[mp],nwave[mp],iexp[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  int nuall,ndall,nu1234,nd5678;
  int nnall,nsall,nn1256,ns3478,nn12,nn56,ns34,ns78;
  int neall,nwall,ne1357,nw2468,ne13,ne57,ne1,ne3,ne5,ne7,nw24,nw68,nw2,nw4,nw6,nw8;
  int nalpha,nalpha2;
  int nl,mth,mmax,ncurrent,nms,ln,lm;
  int iuall[36],ixuall[36],iyuall[36];
  int idall[36],ixdall[36],iydall[36];
  int inall[36],ixnall[36],iynall[36];
  int isall[36],ixsall[36],iysall[36];
  int ieall[36],ixeall[36],iyeall[36];
  int iwall[36],ixwall[36],iywall[36];
  int iu1234[16],ixu1234[16],iyu1234[16];
  int id5678[16],ixd5678[16],iyd5678[16];
  int in1256[16],ixn1256[16],iyn1256[16];
  int is3478[16],ixs3478[16],iys3478[16];
  int in12[16],ixn12[16],iyn12[16];
  int in56[16],ixn56[16],iyn56[16];
  int is34[16],ixs34[16],iys34[16];
  int is78[16],ixs78[16],iys78[16];
  int ie1357[16],ixe1357[16],iye1357[16];
  int iw2468[16],ixw2468[16],iyw2468[16];
  int ie13[16],ixe13[16],iye13[16];
  int ie57[16],ixe57[16],iye57[16];
  int iw24[16],ixw24[16],iyw24[16];
  int iw68[16],ixw68[16],iyw68[16];
  int ie1[16],ixe1[16],iye1[16];
  int ie3[16],ixe3[16],iye3[16];
  int ie5[16],ixe5[16],iye5[16];
  int ie7[16],ixe7[16],iye7[16];
  int iw2[16],ixw2[16],iyw2[16];
  int iw4[16],ixw4[16],iyw4[16];
  int iw6[16],ixw6[16],iyw6[16];
  int iw8[16],ixw8[16],iyw8[16];
  double rb,sgn,halpha,alpha,lscale,rscale;
  double xi[mp],wt[mp];
  std::complex<double> cnm,aaxdd;
  std::complex<double> ztmp1,ztmp2,ztmp,zmul;
  std::complex<double> aax[mpsym],aaxd[mpsym];
  std::complex<double> bbx[mpsym],bbxd[mpsym];
  std::complex<double> exu[162],exd[162],exn[162],exs[162],exe[162],exw[162];
  std::complex<double> fxu[162],fxd[162],fxn[162],fxs[162],fxe[162],fxw[162];
  std::complex<double> fxu1234[162],fxd5678[162];
  std::complex<double> fxuall[162],fxdall[162];
  std::complex<double> fxn12[162],fxn56[162],fxn1256[162],fxnall[162];
  std::complex<double> fxs34[162],fxs78[162],fxs3478[162],fxsall[162];
  std::complex<double> fxe13[162],fxe57[162],fxe1357[162],fxeall[162];
  std::complex<double> fxe1[162],fxe3[162],fxe5[162],fxe7[162];
  std::complex<double> fxw24[162],fxw68[162],fxw2468[162],fxwall[162];
  std::complex<double> fxw2[162],fxw4[162],fxw6[162],fxw8[162];
  std::complex<double> ima(0,1);
  std::complex<double> zeye[200],expl[2000],exmi[2000];
  std::complex<double> ephi[100];
  tic = get_time();

  kl = 1 << lev;
  rb = rd/kl;
  lscale = rd/kl*lambda;

  if( ini != 0 ) {
    for( i=0; i<ini; i++ ) {
      for( j=0; j<nmp; j++ ) {
        ax[i][j] = 0;
      }
    }
  }
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
  precalch(lev,xi,wt,nfour,nwave);
  lexp1 = new std::complex<double>* [nbnet];
  for( i=0; i<nbne; i++ ) lexp1[i] = new std::complex<double> [maxwave];
  lexp2 = new std::complex<double>* [nbnet];
  for( i=0; i<nbne; i++ ) lexp2[i] = new std::complex<double> [maxwave];

// process up and down lists -----------------------------------------------------------------------
  for( ii=0; ii<lbi; ii++ ) {
    ib = ii+nlbi[lev-1];
    iexp1[ib] = 0;
    iexp2[ib] = 0;
    for( jj=0; jj<maxwave; jj++ ) {
      lexp1[ib][jj] = 0;
      lexp2[ib][jj] = 0;
    }
  }

  nbox = -1;
  for( ii=0; ii<lbi; ii++ ) {
    ib = ii+nlbi[lev-1];
    for( i=0; i<nmp; i++ ) {
      bbx[i] = bx[ib][i];
    }

// YMKUDEXP
// YMPOLETOEXP
    ntot = 0;
    for( l=0; l<mp; l++ ) {
      sgn = -1;
      for( m=0; m<nfour[l]; m++ ) {
        ic = ntot+m;
        ztmp1 = 0;
        ztmp2 = 0;
        sgn = -sgn;
        for( n=m; n<mp; n+=2 ) {
          nm = n*(n+1)/2+m;
          ztmp1 += std::complex<double>(m2e[l][m][n])*bbx[nm];
        }
        for( n=m+1; n<mp; n+=2 ) {
          nm = n*(n+1)/2+m;
          ztmp2 += std::complex<double>(m2e[l][m][n])*bbx[nm];
        }
        exu[ic] = ztmp1+ztmp2;
        exd[ic] = sgn*(ztmp1-ztmp2);
      }
      ntot += nfour[l];
    }

// YFTOPHYS
    nftot = 0;
    nptot = 0;
    nexte = 0;
    nexto = 0;
    for( l=0; l<mp; l++ ) {
      for( ival=0; ival<nwave[l]/2; ival++ ) {
        fxu[nptot+ival] = exu[nftot];
        fxd[nptot+ival] = exd[nftot];
        sgn = -2;
        for( mm=1; mm<nfour[l]; mm+=2 ) {
          sgn = -sgn;
          fxu[nptot+ival] += ima*std::complex<double>(sgn*std::real( fexpe[nexte]*exu[nftot+mm] ));
          fxd[nptot+ival] += ima*std::complex<double>(sgn*std::real( fexpe[nexte]*exd[nftot+mm] ));
          nexte++;
        }
        sgn = 2;
        for( mm=2; mm<nfour[l]; mm+=2 ) {
          sgn = -sgn;
          fxu[nptot+ival] += sgn*std::real( fexpo[nexto]*exu[nftot+mm] );;
          fxd[nptot+ival] += sgn*std::real( fexpo[nexto]*exd[nftot+mm] );;
          nexto++;
        }
      }
      nftot += nfour[l];
      nptot += nwave[l]/2;
    }

// YMKUDEXP part2
    islast = 1;
    if( nbox != nfi[ib]/8 ) {
      nbox = nfi[ib]/8;
      iexp[0]=0;
      iexp[1]=0;
      iexp[2]=0;
      iexp[3]=0;
      for( jj=0; jj<maxwave; jj++ ) {
        fxu1234[jj] = 0;
        fxd5678[jj] = 0;
        fxuall[jj] = 0;
        fxdall[jj] = 0;
      }
    }

    if( nfi[ib] % 8 == 0 ) {
      for( jj=0; jj<maxwave; jj++ ) {
        fxu1234[jj] += fxu[jj];
        fxdall[jj] += fxd[jj];
      }
      iexp[1]++;
      iexp[2]++;
    } else if( nfi[ib] % 8 == 1 ) {
      for( jj=0; jj<maxwave; jj++ ) {
        fxu[jj] *= std::conj(xs[jj][0]);
        fxu1234[jj] += fxu[jj];
        fxd[jj] *= xs[jj][0];
        fxdall[jj] += fxd[jj];
      }
      iexp[1]++;
      iexp[2]++;
    } else if( nfi[ib] % 8 == 2 ) {
      for( jj=0; jj<maxwave; jj++ ) {
        fxu[jj] *= std::conj(ys[jj][0]);
        fxu1234[jj] += fxu[jj];
        fxd[jj] *= ys[jj][0];
        fxdall[jj] += fxd[jj];
      }
      iexp[1]++;
      iexp[2]++;
    } else if( nfi[ib] % 8 == 3 ) {
      for( jj=0; jj<maxwave; jj++ ) {
        ztmp = xs[jj][0]*ys[jj][0];
        fxu[jj] *= std::conj(ztmp);
        fxu1234[jj] += fxu[jj];
        fxd[jj] *= ztmp;
        fxdall[jj] += fxd[jj];
      }
      iexp[1]++;
      iexp[2]++;
    } else if( nfi[ib] % 8 == 4 ) {
      for( jj=0; jj<maxwave; jj++ ) {
        fxu[jj] /= zs[jj][0];
        fxuall[jj] += fxu[jj];
        fxd[jj] *= zs[jj][0];
        fxd5678[jj] += fxd[jj];
      }
      iexp[0]++;
      iexp[3]++;
    } else if( nfi[ib] % 8 == 5 ) {
      for( jj=0; jj<maxwave; jj++ ) {
        ztmp = xs[jj][0]*zs[jj][0];
        fxu[jj] /= ztmp;
        fxuall[jj] += fxu[jj];
        fxd[jj] *= ztmp;
        fxd5678[jj] += fxd[jj];
      }
      iexp[0]++;
      iexp[3]++;
    } else if( nfi[ib] % 8 == 6 ) {
      for( jj=0; jj<maxwave; jj++ ) {
        ztmp = ys[jj][0]*zs[jj][0];
        fxu[jj] /= ztmp;
        fxuall[jj] += fxu[jj];
        fxd[jj] *= ztmp;
        fxd5678[jj] += fxd[jj];
      }
      iexp[0]++;
      iexp[3]++;
    } else if( nfi[ib] % 8 == 7 ) {
      for( jj=0; jj<maxwave; jj++ ) {
        ztmp = xs[jj][0]*ys[jj][0]*zs[jj][0];
        fxu[jj] /= ztmp;
        fxuall[jj] += fxu[jj];
        fxd[jj] *= ztmp;
        fxd5678[jj] += fxd[jj];
      }
      iexp[0]++;
      iexp[3]++;
    }
    if( ib == lbi-1 ) {
      for( jj=0; jj<maxwave; jj++ ) {
        fxuall[jj] += fxu1234[jj];
        fxdall[jj] += fxd5678[jj];
      }
      islast = 0;
    } else if( nfi[ib+1]/8 != nbox ) {
      for( jj=0; jj<maxwave; jj++ ) {
        fxuall[jj] += fxu1234[jj];
        fxdall[jj] += fxd5678[jj];
      }
      islast = 0;
    }

    if( islast == 0 ) {

// MKUPLIST
      nuall = 0;
      nu1234 = 0;
      boxc(nfi[ib]/8,3,nc);
      ixp = nc[0]+1;
      iyp = nc[1]+1;
      izp = nc[2]+1;
      jzp = izp+1;
      for( jyp=iyp-1; jyp<=iyp+1; jyp++ ) {
        for( jxp=ixp-1; jxp<=ixp+1; jxp++ ) {
          for( jz=std::max(2*jzp-2,jzmin); jz<=std::min(2*jzp-1,jzmax); jz++ ) {
            for( jy=std::max(2*jyp-2,jymin); jy<=std::min(2*jyp-1,jymax); jy++ ) {
              for( jx=std::max(2*jxp-2,jxmin); jx<=std::min(2*jxp-1,jxmax); jx++ ) {
                ix = jx-2*(ixp-1);
                iy = jy-2*(iyp-1);
                iz = jz-2*(izp-1);
                if( iz > 2.5 ) {
                  nc[0] = jx;
                  nc[1] = jy;
                  nc[2] = jz;
                  boxn1(nc,je,lev);
                  jj = nej[je];
                  if( jj != -1 ) {
                    iuall[nuall] = jj;
                    ixuall[nuall] = ix;
                    iyuall[nuall] = iy;
                    nuall++;
                  }
                } else if( ix != -2 && ix != 3 && iy != -2 && iy != 3 ) {
                  nc[0] = jx;
                  nc[1] = jy;
                  nc[2] = jz;
                  boxn1(nc,je,lev);
                  jj = nej[je];
                  if( jj != -1 ) {
                    iu1234[nu1234] = jj;
                    ixu1234[nu1234] = ix;
                    iyu1234[nu1234] = iy;
                    nu1234++;
                  }
                }
              }
            }
          }
        }
      }

// YPROCESSUP
      if( iexp[0] > 0 ) {
        for( i=0; i<nuall; i++ ) {
          iexp1[iuall[i]]++;
          for( jj=0; jj<maxwave; jj++ ) {
            zmul = zs[jj][2];
            if( ixuall[i] > 0 ) zmul *= xs[jj][ixuall[i]-1];
            if( ixuall[i] < 0 ) zmul *= std::conj(xs[jj][-ixuall[i]-1]);
            if( iyuall[i] > 0 ) zmul *= ys[jj][iyuall[i]-1];
            if( iyuall[i] < 0 ) zmul *= std::conj(ys[jj][-iyuall[i]-1]);
            lexp1[iuall[i]][jj] += fxuall[jj]*zmul;
          }
        }
      }
      if( iexp[1] > 0 ) {
        for( i=0; i<nu1234; i++ ) {
          iexp1[iu1234[i]]++;
          for( jj=0; jj<maxwave; jj++ ) {
            zmul = zs[jj][1];
            if( ixu1234[i] > 0 ) zmul *= xs[jj][ixu1234[i]-1];
            if( ixu1234[i] < 0 ) zmul *= std::conj(xs[jj][-ixu1234[i]-1]);
            if( iyu1234[i] > 0 ) zmul *= ys[jj][iyu1234[i]-1];
            if( iyu1234[i] < 0 ) zmul *= std::conj(ys[jj][-iyu1234[i]-1]);
            lexp1[iu1234[i]][jj] += fxu1234[jj]*zmul;
          }
        }
      }

// MKDNLIST
      ndall = 0;
      nd5678 = 0;
      boxc(nfi[ib]/8,3,nc);
      ixp = nc[0]+1;
      iyp = nc[1]+1;
      izp = nc[2]+1;
      jzp = izp-1;
      for( jyp=iyp-1; jyp<=iyp+1; jyp++ ) {
        for( jxp=ixp-1; jxp<=ixp+1; jxp++ ) {
          for( jz=std::max(2*jzp-2,jzmin); jz<=std::min(2*jzp-1,jzmax); jz++ ) {
            for( jy=std::max(2*jyp-2,jymin); jy<=std::min(2*jyp-1,jymax); jy++ ) {
              for( jx=std::max(2*jxp-2,jxmin); jx<=std::min(2*jxp-1,jxmax); jx++ ) {
                ix = jx-2*(ixp-1);
                iy = jy-2*(iyp-1);
                iz = jz-2*(izp-1);
                if( iz <= -1.5 ) {
                  nc[0] = jx;
                  nc[1] = jy;
                  nc[2] = jz;
                  boxn1(nc,je,lev);
                  jj = nej[je];
                  if( jj != -1 ) {
                    idall[ndall] = jj;
                    ixdall[ndall] = ix;
                    iydall[ndall] = iy;
                    ndall++;
                  }
                } else if( ix != -2 && ix != 3 && iy != -2 && iy != 3 ){
                  nc[0] = jx;
                  nc[1] = jy;
                  nc[2] = jz;
                  boxn1(nc,je,lev);
                  jj = nej[je];
                  if( jj != -1 ) {
                    id5678[nd5678] = jj;
                    ixd5678[nd5678] = ix;
                    iyd5678[nd5678] = iy;
                    nd5678++;
                  }
                }
              }
            }
          }
        }
      }

// YPROCESSDN
      if( iexp[2] > 0 ) {
        for( i=0; i<ndall; i++ ) {
          iexp2[idall[i]]++;
          for( jj=0; jj<maxwave; jj++ ) {
            zmul = zs[jj][1];
            if( ixdall[i] > 0 ) zmul *= std::conj(xs[jj][ixdall[i]-1]);
            if( ixdall[i] < 0 ) zmul *= xs[jj][-ixdall[i]-1];
            if( iydall[i] > 0 ) zmul *= std::conj(ys[jj][iydall[i]-1]);
            if( iydall[i] < 0 ) zmul *= ys[jj][-iydall[i]-1];
            lexp2[idall[i]][jj] += fxdall[jj]*zmul;
          }
        }
      }
      if( iexp[3] > 0 ) {
        for( i=0; i<nd5678; i++ ) {
          iexp2[id5678[i]]++;
          for( jj=0; jj<maxwave; jj++ ) {
            zmul = zs[jj][0];
            if( ixd5678[i] > 0 ) zmul *= std::conj(xs[jj][ixd5678[i]-1]);
            if( ixd5678[i] < 0 ) zmul *= xs[jj][-ixd5678[i]-1];
            if( iyd5678[i] > 0 ) zmul *= std::conj(ys[jj][iyd5678[i]-1]);
            if( iyd5678[i] < 0 ) zmul *= ys[jj][-iyd5678[i]-1];
            lexp2[id5678[i]][jj] += fxd5678[jj]*zmul;
          }
        }
      }
    }
  }

  for( ii=0; ii<lbi; ii++ ) {
    if( iexp1[ii] > 0 ) {
// YPHYSTOF
      nftot = 0;
      nptot = 0;
      next = 0;
      for( l=0; l<mp; l++ ) {
        nalpha = nwave[l];
        nalpha2 = nalpha/2;
        halpha = 2*pi/nalpha;
        exu[nftot] = 0;
        for( ival=0; ival<nalpha2; ival++ ) {
          exu[nftot] += 2.0*std::real(lexp1[ii][nptot+ival]);
        }
        exu[nftot] /= nalpha;
        for( mm=2; mm<nfour[l]; mm+=2 ) {
          exu[nftot+mm] = 0;
          for( ival=0; ival<nalpha2; ival++ ) {
            exu[nftot+mm] += fexpb[next]*2.0*std::real(lexp1[ii][nptot+ival]);
            next++;
          }
          exu[nftot+mm] /= nalpha;
        }
        for( mm=1; mm<nfour[l]; mm+=2 ) {
          exu[nftot+mm] = 0;
          for( ival=0; ival<nalpha2; ival++ ) {
            exu[nftot+mm] += fexpb[next]*2.0*ima*std::complex<double>(std::imag(lexp1[ii][nptot+ival]));
            next++;
          }
          exu[nftot+mm] /= nalpha;
        }
        nftot += nfour[l];
        nptot += nwave[l]/2;
      }
    }
    if( iexp2[ii] > 0 ) {
// YPHYSTOF
      nftot = 0;
      nptot = 0;
      next = 0;
      for( l=0; l<mp; l++ ) {
        nalpha = nwave[l];
        nalpha2 = nalpha/2;
        halpha = 2*pi/nalpha;
        exd[nftot] = 0;
        for( ival=0; ival<nalpha2; ival++ ) {
          exd[nftot] += 2.0*std::real(lexp2[ii][nptot+ival]);
        }
        exd[nftot] /= nalpha;
        for( mm=2; mm<nfour[l]; mm+=2 ) {
          exd[nftot+mm] = 0;
          for( ival=0; ival<nalpha2; ival++ ) {
            exd[nftot+mm] += fexpb[next]*2.0*std::real(lexp2[ii][nptot+ival]);
            next++;
          }
          exd[nftot+mm] /= nalpha;
        }
        for( mm=1; mm<nfour[l]; mm+=2 ) {
          exd[nftot+mm] = 0;
          for( ival=0; ival<nalpha2; ival++ ) {
            exd[nftot+mm] += fexpb[next]*2.0*ima*std::complex<double>(std::imag(lexp2[ii][nptot+ival]));
            next++;
          }
          exd[nftot+mm] /= nalpha;
        }
        nftot += nfour[l];
        nptot += nwave[l]/2;
      }
    }
    if( iexp1[ii] > 0 || iexp2[ii] > 0 ) {
// YEXPTOLOCAL
      zeye[0] = 1;
      for( i=1; i<nmax; i++ ) {
        zeye[i] = zeye[i-1]*ima;
      }
      for( nm=0; nm<mp; nm++ ) {
        for( mth=0; mth<mp; mth++ ) {
          nms = nm*(nm+1)/2+mth;
          aax[nms] = 0;
        }
      }
      if( iexp1[ii] <= 0 ) {
        for( nm=0; nm<maxfour; nm++ ) {
          expl[nm] = exd[nm];
          exmi[nm] = exd[nm];
        }
      } else if( iexp2[ii] <= 0 ) {
        for( nm=0; nm<maxfour; nm++ ) {
          expl[nm] = exu[nm];
          exmi[nm] = -exu[nm];
        }
      } else {
        for( nm=0; nm<maxfour; nm++ ) {
          expl[nm] = exd[nm] + exu[nm];
          exmi[nm] = exd[nm] - exu[nm];
        }
      }
      ntot = 0;
      for( nl=0; nl<mp; nl++ ) {
        for( nm=0; nm<mp; nm+=2 ) {
          mmax = nfour[nl]-1;
          if( mmax > nm ) mmax = nm;
          for( mth=0; mth<=mmax; mth++ ) {
            nms = nm*(nm+1)/2+mth;
            ncurrent = ntot+mth;
            aax[nms] += m2e[nl][mth][nm]*wt[nl]*expl[ncurrent];
          }
        }
        for( nm=1; nm<mp; nm+=2 ) {
          mmax = nfour[nl]-1;
          if( mmax > nm ) mmax = nm;
          for( mth=0; mth<=mmax; mth++ ) {
            nms = nm*(nm+1)/2+mth;
            ncurrent = ntot+mth;
            aax[nms] += m2e[nl][mth][nm]*wt[nl]*exmi[ncurrent];
          }
        }
        ntot += nfour[nl];
      }
      rscale = pi/lscale*0.5;
      for( nm=0; nm<mp; nm++ ) {
        for( mth=0; mth<=nm; mth++ ) {
          nms = nm*(nm+1)/2+mth;
          aax[nms] *= zeye[mth]*rscale*double(ytop[nm][mth]);
          ax[ii][nms] += aax[nms];
//          if(fabs(std::real(expl[0])-113.47)<0.01) std::cout << nm << " " << mth << " " << aax[nms] << std::endl;
        }
      }
    }
  }

// process north and south lists -------------------------------------------------------------------
  for( ii=0; ii<lbi; ii++ ) {
    ib = ii+nlbi[lev-1];
    iexp1[ib] = 0;
    iexp2[ib] = 0;
    for( jj=0; jj<maxwave; jj++ ) {
      lexp1[ib][jj] = 0;
      lexp2[ib][jj] = 0;
    }
  }

  nbox = -1;
  for( ii=0; ii<lbi; ii++ ) {
    ib = ii+nlbi[lev-1];
    for( i=0; i<nmp; i++ ) {
      bbxd[i] = bx[ib][i];
    }

// ROTZTOY
    ephi[0] = 1;
    for( m=1; m<mp+1; m++ ) {
      ephi[m] = ephi[m-1]*std::conj(ima);
    }
    for( l=0; l<mp; l++ ) {
      for( m=0; m<=l; m++ ) {
        ln = l*(l+1)/2;
        lm = l*(l+1)/2+m;
        bbx[lm] = ephi[0]*bbxd[ln]*rdmpi2[l][0][mp+m];
        for( n=1; n<=l; n++ ) {
          ln = l*(l+1)/2+n;
          bbx[lm] += ephi[n]*bbxd[ln]*rdmpi2[l][n][mp+m]+std::conj(ephi[n]*bbxd[ln])*rdmpi2[l][n][mp-m];
        }
      }
    }

// YMKNSEXP
// YMPOLETOEXP
    ntot = 0;
    for( l=0; l<mp; l++ ) {
      sgn = -1;
      for( m=0; m<nfour[l]; m++ ) {
        ic = ntot+m;
        ztmp1 = 0;
        ztmp2 = 0;
        sgn = -sgn;
        for( n=m; n<mp; n+=2 ) {
          nm = n*(n+1)/2+m;
          ztmp1 += std::complex<double>(m2e[l][m][n])*bbx[nm];
        }
        for( n=m+1; n<mp; n+=2 ) {
          nm = n*(n+1)/2+m;
          ztmp2 += std::complex<double>(m2e[l][m][n])*bbx[nm];
        }
        exn[ic] = ztmp1+ztmp2;
        exs[ic] = sgn*(ztmp1-ztmp2);
      }
      ntot += nfour[l];
    }

// YFTOPHYS
    nftot = 0;
    nptot = 0;
    nexte = 0;
    nexto = 0;
    for( l=0; l<mp; l++ ) {
      for( ival=0; ival<nwave[l]/2; ival++ ) {
        fxn[nptot+ival] = exn[nftot];
        fxs[nptot+ival] = exs[nftot];
        sgn = -2;
        for( mm=1; mm<nfour[l]; mm+=2 ) {
          sgn = -sgn;
          fxn[nptot+ival] += ima*std::complex<double>(sgn*std::real( fexpe[nexte]*exn[nftot+mm] ));
          fxs[nptot+ival] += ima*std::complex<double>(sgn*std::real( fexpe[nexte]*exs[nftot+mm] ));
          nexte++;
        }
        sgn = 2;
        for( mm=2; mm<nfour[l]; mm+=2 ) {
          sgn = -sgn;
          fxn[nptot+ival] += sgn*std::real( fexpo[nexto]*exn[nftot+mm] );;
          fxs[nptot+ival] += sgn*std::real( fexpo[nexto]*exs[nftot+mm] );;
          nexto++;
        }
      }
      nftot += nfour[l];
      nptot += nwave[l]/2;
    }

// YMKNSEXP part2
    islast = 1;
    if( nbox != nfi[ib]/8 ) {
      nbox = nfi[ib]/8;
      iexp[0]=0;
      iexp[1]=0;
      iexp[2]=0;
      iexp[3]=0;
      iexp[4]=0;
      iexp[5]=0;
      iexp[6]=0;
      iexp[7]=0;
      for( jj=0; jj<maxwave; jj++ ) {
        fxnall[jj] = 0;
        fxn1256[jj] = 0;
        fxn12[jj] = 0;
        fxn56[jj] = 0;
        fxsall[jj] = 0;
        fxs3478[jj] = 0;
        fxs34[jj] = 0;
        fxs78[jj] = 0;
      }
    }

    if( nfi[ib] % 8 == 0 ) {
      for( jj=0; jj<maxwave; jj++ ) {
        fxn12[jj] += fxn[jj];
        fxsall[jj] += fxs[jj];
      }
      iexp[2]++;
      iexp[4]++;
    } else if( nfi[ib] % 8 == 1 ) {
      for( jj=0; jj<maxwave; jj++ ) {
        fxn[jj] *= std::conj(ys[jj][0]);
        fxn12[jj] += fxn[jj];
        fxs[jj] *= ys[jj][0];
        fxsall[jj] += fxs[jj];
      }
      iexp[2]++;
      iexp[4]++;
    } else if( nfi[ib] % 8 == 2 ) {
      for( jj=0; jj<maxwave; jj++ ) {
        fxn[jj] /= zs[jj][0];
        fxnall[jj] += fxn[jj];
        fxs[jj] *= zs[jj][0];
        fxs34[jj] += fxs[jj];
      }
      iexp[0]++;
      iexp[6]++;
    } else if( nfi[ib] % 8 == 3 ) {
      for( jj=0; jj<maxwave; jj++ ) {
        ztmp = ys[jj][0]*zs[jj][0];
        fxn[jj] /= ztmp;
        fxnall[jj] += fxn[jj];
        fxs[jj] *= ztmp;
        fxs34[jj] += fxs[jj];
      }
      iexp[0]++;
      iexp[6]++;
    } else if( nfi[ib] % 8 == 4 ) {
      for( jj=0; jj<maxwave; jj++ ) {
        fxn[jj] *= std::conj(xs[jj][0]);
        fxn56[jj] += fxn[jj];
        fxs[jj] *= xs[jj][0];
        fxsall[jj] += fxs[jj];
      }
      iexp[3]++;
      iexp[4]++;
    } else if( nfi[ib] % 8 == 5 ) {
      for( jj=0; jj<maxwave; jj++ ) {
        ztmp = ys[jj][0]*xs[jj][0];
        fxn[jj] *= std::conj(ztmp);
        fxn56[jj] += fxn[jj];
        fxs[jj] *= ztmp;
        fxsall[jj] += fxs[jj];
      }
      iexp[3]++;
      iexp[4]++;
    } else if( nfi[ib] % 8 == 6 ) {
      for( jj=0; jj<maxwave; jj++ ) {
        ztmp = zs[jj][0]*xs[jj][0];
        fxn[jj] /= ztmp;
        fxnall[jj] += fxn[jj];
        fxs[jj] *= ztmp;
        fxs78[jj] += fxs[jj];
      }
      iexp[0]++;
      iexp[7]++;
    } else if( nfi[ib] % 8 == 7 ) {
      for( jj=0; jj<maxwave; jj++ ) {
        ztmp = xs[jj][0]*ys[jj][0]*zs[jj][0];
        fxn[jj] /= ztmp;
        fxnall[jj] += fxn[jj];
        fxs[jj] *= ztmp;
        fxs78[jj] += fxs[jj];
      }
      iexp[0]++;
      iexp[7]++;
    }
    if( ib == lbi-1 ) {
      for( jj=0; jj<maxwave; jj++ ) {
        fxn1256[jj] = fxn56[jj]+fxn12[jj];
        fxs3478[jj] = fxs78[jj]+fxs34[jj];
        fxnall[jj] += fxn1256[jj];
        fxsall[jj] += fxs3478[jj];
      }
      islast = 0;
      iexp[1]++;
      iexp[4]++;
      iexp[5]++;
    } else if ( nfi[ib+1]/8 != nbox ) {
      for( jj=0; jj<maxwave; jj++ ) {
        fxn1256[jj] = fxn56[jj]+fxn12[jj];
        fxs3478[jj] = fxs78[jj]+fxs34[jj];
        fxnall[jj] += fxn1256[jj];
        fxsall[jj] += fxs3478[jj];
      }
      islast = 0;
      iexp[1]++;
      iexp[4]++;
      iexp[5]++;
    }

    if( islast == 0 ) {

// MKNOLIST
      nnall = 0;
      nn1256 = 0;
      nn12 = 0;
      nn56 = 0;
      boxc(nfi[ib]/8,3,nc);
      ixp = nc[0]+1;
      iyp = nc[1]+1;
      izp = nc[2]+1;
      jyp = iyp+1;
      for( jzp=izp-1; jzp<=izp+1; jzp++ ) {
        for( jxp=ixp-1; jxp<=ixp+1; jxp++ ) {
          for( jz=std::max(2*jzp-2,jzmin); jz<=std::min(2*jzp-1,jzmax); jz++ ) {
            for( jy=std::max(2*jyp-2,jymin); jy<=std::min(2*jyp-1,jymax); jy++ ) {
              for( jx=std::max(2*jxp-2,jxmin); jx<=std::min(2*jxp-1,jxmax); jx++ ) {
                ix = jx-2*(ixp-1);
                iy = jy-2*(iyp-1);
                iz = jz-2*(izp-1);
                if( iy > 2.5 ) {
                  if( iz != -2 && iz != 3 ) {
                    nc[0] = jx;
                    nc[1] = jy;
                    nc[2] = jz;
                    boxn1(nc,je,lev);
                    jj = nej[je];
                    if( jj != -1 ) {
                      inall[nnall] = jj;
                      ixnall[nnall] = iz;
                      iynall[nnall] = ix;
                      nnall++;
                    }
                  }
                } else if( ( iz == 0 || iz == 1 ) && ix >= -1 && ix <= 2 ) {
                  nc[0] = jx;
                  nc[1] = jy;
                  nc[2] = jz;
                  boxn1(nc,je,lev);
                  jj = nej[je];
                  if( jj != -1 ) {
                    in1256[nn1256] = jj;
                    ixn1256[nn1256] = iz;
                    iyn1256[nn1256] = ix;
                    nn1256++;
                  }
                } else if( iz == -1 && ( ix >= -1 && ix <= 2) ) {
                  nc[0] = jx;
                  nc[1] = jy;
                  nc[2] = jz;
                  boxn1(nc,je,lev);
                  jj = nej[je];
                  if( jj != -1 ) {
                    in12[nn12] = jj;
                    ixn12[nn12] = iz;
                    iyn12[nn12] = ix;
                    nn12++;
                  }
                } else if( iz == 2 && ( ix >= -1 && ix <= 2 ) ) {
                  nc[0] = jx;
                  nc[1] = jy;
                  nc[2] = jz;
                  boxn1(nc,je,lev);
                  jj = nej[je];
                  if( jj != -1 ) {
                    in56[nn56] = jj;
                    ixn56[nn56] = iz;
                    iyn56[nn56] = ix;
                    nn56++;
                  }
                }
              }
            }
          }
        }
      }

// YPROCESSNO
      if( iexp[0] > 0 ) {
        for( i=0; i<nnall; i++ ) {
          iexp1[inall[i]]++;
          for( jj=0; jj<maxwave; jj++ ) {
            zmul = zs[jj][2];
            if( ixnall[i] > 0 ) zmul *= xs[jj][ixnall[i]-1];
            if( ixnall[i] < 0 ) zmul *= std::conj(xs[jj][-ixnall[i]-1]);
            if( iynall[i] > 0 ) zmul *= ys[jj][iynall[i]-1];
            if( iynall[i] < 0 ) zmul *= std::conj(ys[jj][-iynall[i]-1]);
            lexp1[inall[i]][jj] += fxnall[jj]*zmul;
          }
        }
      }
      if( iexp[1] > 0 ) {
        for( i=0; i<nn1256; i++ ) {
          iexp1[in1256[i]]++;
          for( jj=0; jj<maxwave; jj++ ) {
            zmul = zs[jj][1];
            if( ixn1256[i] > 0 ) zmul *= xs[jj][ixn1256[i]-1];
            if( ixn1256[i] < 0 ) zmul *= std::conj(xs[jj][-ixn1256[i]-1]);
            if( iyn1256[i] > 0 ) zmul *= ys[jj][iyn1256[i]-1];
            if( iyn1256[i] < 0 ) zmul *= std::conj(ys[jj][-iyn1256[i]-1]);
            lexp1[in1256[i]][jj] += fxn1256[jj]*zmul;
          }
        }
      }
      if( iexp[2] > 0 ) {
        for( i=0; i<nn12; i++ ) {
          iexp1[in12[i]]++;
          for( jj=0; jj<maxwave; jj++ ) {
            zmul = zs[jj][1];
            if( ixn12[i] > 0 ) zmul *= xs[jj][ixn12[i]-1];
            if( ixn12[i] < 0 ) zmul *= std::conj(xs[jj][-ixn12[i]-1]);
            if( iyn12[i] > 0 ) zmul *= ys[jj][iyn12[i]-1];
            if( iyn12[i] < 0 ) zmul *= std::conj(ys[jj][-iyn12[i]-1]);
            lexp1[in12[i]][jj] += fxn12[jj]*zmul;
          }
        }
      }
      if( iexp[3] > 0 ) {
        for( i=0; i<nn56; i++ ) {
          iexp1[in56[i]]++;
          for( jj=0; jj<maxwave; jj++ ) {
            zmul = zs[jj][1];
            if( ixn56[i] > 0 ) zmul *= xs[jj][ixn56[i]-1];
            if( ixn56[i] < 0 ) zmul *= std::conj(xs[jj][-ixn56[i]-1]);
            if( iyn56[i] > 0 ) zmul *= ys[jj][iyn56[i]-1];
            if( iyn56[i] < 0 ) zmul *= std::conj(ys[jj][-iyn56[i]-1]);
            lexp1[in56[i]][jj] += fxn56[jj]*zmul;
          }
        }
      }

// MKSOLIST
      nsall = 0;
      ns3478 = 0;
      ns34 = 0;
      ns78 = 0;
      boxc(nfi[ib]/8,3,nc);
      ixp = nc[0]+1;
      iyp = nc[1]+1;
      izp = nc[2]+1;
      jyp = iyp-1;
      for( jzp=izp-1; jzp<=izp+1; jzp++ ) {
        for( jxp=ixp-1; jxp<=ixp+1; jxp++ ) {
          for( jz=std::max(2*jzp-2,jzmin); jz<=std::min(2*jzp-1,jzmax); jz++ ) {
            for( jy=std::max(2*jyp-2,jymin); jy<=std::min(2*jyp-1,jymax); jy++ ) {
              for( jx=std::max(2*jxp-2,jxmin); jx<=std::min(2*jxp-1,jxmax); jx++ ) {
                ix = jx-2*(ixp-1);
                iy = jy-2*(iyp-1);
                iz = jz-2*(izp-1);
                if( iy <= -1.5 ) {
                  if( iz != -2 && iz != 3 ) {
                    nc[0] = jx;
                    nc[1] = jy;
                    nc[2] = jz;
                    boxn1(nc,je,lev);
                    jj = nej[je];
                    if( jj != -1 ) {
                      isall[nsall] = jj;
                      ixsall[nsall] = iz;
                      iysall[nsall] = ix;
                      nsall++;
                    }
                  }
                } else if( ( iz == 0 || iz == 1 ) && ix >= -1 && ix <= 2 ){
                  nc[0] = jx;
                  nc[1] = jy;
                  nc[2] = jz;
                  boxn1(nc,je,lev);
                  jj = nej[je];
                  if( jj != -1 ) {
                    is3478[ns3478] = jj;
                    ixs3478[ns3478] = iz;
                    iys3478[ns3478] = ix;
                    ns3478++;
                  }
                } else if( iz == -1 && ( ix >= -1 && ix <= 2 ) ) {
                  nc[0] = jx;
                  nc[1] = jy;
                  nc[2] = jz;
                  boxn1(nc,je,lev);
                  jj = nej[je];
                  if( jj != -1 ) {
                    is34[ns34] = jj;
                    ixs34[ns34] = iz;
                    iys34[ns34] = ix;
                    ns34++;
                  }
                } else if( iz == 2 && ( ix >= -1 && ix <= 2 ) ) {
                  nc[0] = jx;
                  nc[1] = jy;
                  nc[2] = jz;
                  boxn1(nc,je,lev);
                  jj = nej[je];
                  if( jj != -1 ) {
                    is78[ns78] = jj;
                    ixs78[ns78] = iz;
                    iys78[ns78] = ix;
                    ns78++;
                  }
                }
              }
            }
          }
        }
      }

// YPROCESSSO
      if( iexp[4] > 0 ) {
        for( i=0; i<nsall; i++ ) {
          iexp2[isall[i]]++;
          for( jj=0; jj<maxwave; jj++ ) {
            zmul = zs[jj][1];
            if( ixsall[i] > 0 ) zmul *= std::conj(xs[jj][ixsall[i]-1]);
            if( ixsall[i] < 0 ) zmul *= xs[jj][-ixsall[i]-1];
            if( iysall[i] > 0 ) zmul *= std::conj(ys[jj][iysall[i]-1]);
            if( iysall[i] < 0 ) zmul *= ys[jj][-iysall[i]-1];
            lexp2[isall[i]][jj] += fxsall[jj]*zmul;
          }
        }
      }
      if( iexp[5] > 0 ) {
        for( i=0; i<ns3478; i++ ) {
          iexp2[is3478[i]]++;
          for( jj=0; jj<maxwave; jj++ ) {
            zmul = zs[jj][0];
            if( ixs3478[i] > 0 ) zmul *= std::conj(xs[jj][ixs3478[i]-1]);
            if( ixs3478[i] < 0 ) zmul *= xs[jj][-ixs3478[i]-1];
            if( iys3478[i] > 0 ) zmul *= std::conj(ys[jj][iys3478[i]-1]);
            if( iys3478[i] < 0 ) zmul *= ys[jj][-iys3478[i]-1];
            lexp2[is3478[i]][jj] += fxs3478[jj]*zmul;
          }
        }
      }
      if( iexp[6] > 0 ) {
        for( i=0; i<ns34; i++ ) {
          iexp2[is34[i]]++;
          for( jj=0; jj<maxwave; jj++ ) {
            zmul = zs[jj][0];
            if( ixs34[i] > 0 ) zmul *= std::conj(xs[jj][ixs34[i]-1]);
            if( ixs34[i] < 0 ) zmul *= xs[jj][-ixs34[i]-1];
            if( iys34[i] > 0 ) zmul *= std::conj(ys[jj][iys34[i]-1]);
            if( iys34[i] < 0 ) zmul *= ys[jj][-iys34[i]-1];
            lexp2[is34[i]][jj] += fxs34[jj]*zmul;
          }
        }
      }
      if( iexp[7] > 0 ) {
        for( i=0; i<ns78; i++ ) {
          iexp2[is78[i]]++;
          for( jj=0; jj<maxwave; jj++ ) {
            zmul = zs[jj][0];
            if( ixs78[i] > 0 ) zmul *= std::conj(xs[jj][ixs78[i]-1]);
            if( ixs78[i] < 0 ) zmul *= xs[jj][-ixs78[i]-1];
            if( iys78[i] > 0 ) zmul *= std::conj(ys[jj][iys78[i]-1]);
            if( iys78[i] < 0 ) zmul *= ys[jj][-iys78[i]-1];
            lexp2[is78[i]][jj] += fxs78[jj]*zmul;
          }
        }
      }
    }
  }

  for( ii=0; ii<lbi; ii++ ) {
    if( iexp1[ii] > 0 ) {
// YPHYSTOF
      nftot = 0;
      nptot = 0;
      next = 0;
      for( l=0; l<mp; l++ ) {
        nalpha = nwave[l];
        nalpha2 = nalpha/2;
        halpha = 2*pi/nalpha;
        exn[nftot] = 0;
        for( ival=0; ival<nalpha2; ival++ ) {
          exn[nftot] += 2.0*std::real(lexp1[ii][nptot+ival]);
        }
        exn[nftot] /= nalpha;
        for( mm=2; mm<nfour[l]; mm+=2 ) {
          exn[nftot+mm] = 0;
          for( ival=0; ival<nalpha2; ival++ ) {
            exn[nftot+mm] += fexpb[next]*2.0*std::real(lexp1[ii][nptot+ival]);
            next++;
          }
          exn[nftot+mm] /= nalpha;
        }
        for( mm=1; mm<nfour[l]; mm+=2 ) {
          exn[nftot+mm] = 0;
          for( ival=0; ival<nalpha2; ival++ ) {
            exn[nftot+mm] += fexpb[next]*2.0*ima*std::complex<double>(std::imag(lexp1[ii][nptot+ival]));
            next++;
          }
          exn[nftot+mm] /= nalpha;
        }
        nftot += nfour[l];
        nptot += nwave[l]/2;
      }
    }
    if( iexp2[ii] > 0 ) {
// YPHYSTOF
      nftot = 0;
      nptot = 0;
      next = 0;
      for( l=0; l<mp; l++ ) {
        nalpha = nwave[l];
        nalpha2 = nalpha/2;
        halpha = 2*pi/nalpha;
        exs[nftot] = 0;
        for( ival=0; ival<nalpha2; ival++ ) {
          exs[nftot] += 2.0*std::real(lexp2[ii][nptot+ival]);
        }
        exs[nftot] /= nalpha;
        for( mm=2; mm<nfour[l]; mm+=2 ) {
          exs[nftot+mm] = 0;
          for( ival=0; ival<nalpha2; ival++ ) {
            exs[nftot+mm] += fexpb[next]*2.0*std::real(lexp2[ii][nptot+ival]);
            next++;
          }
          exs[nftot+mm] /= nalpha;
        }
        for( mm=1; mm<nfour[l]; mm+=2 ) {
          exs[nftot+mm] = 0;
          for( ival=0; ival<nalpha2; ival++ ) {
            exs[nftot+mm] += fexpb[next]*2.0*ima*std::complex<double>(std::imag(lexp2[ii][nptot+ival]));
            next++;
          }
          exs[nftot+mm] /= nalpha;
        }
        nftot += nfour[l];
        nptot += nwave[l]/2;
      }
    }
    if( iexp1[ii] > 0 || iexp2[ii] > 0 ) {
// YEXPTOLOCAL
      zeye[0] = 1;
      for( i=1; i<nmax; i++ ) {
        zeye[i] = zeye[i-1]*ima;
      }
      for( nm=0; nm<mp; nm++ ) {
        for( mth=0; mth<mp; mth++ ) {
          nms = nm*(nm+1)/2+mth;
          aaxd[nms] = 0;
        }
      }
      if( iexp1[ii] <= 0 ) {
        for( nm=0; nm<maxfour; nm++ ) {
          expl[nm] = exs[nm];
          exmi[nm] = exs[nm];
        }
      } else if( iexp2[ii] <= 0 ) {
        for( nm=0; nm<maxfour; nm++ ) {
          expl[nm] = exn[nm];
          exmi[nm] = -exn[nm];
        }
      } else {
        for( nm=0; nm<maxfour; nm++ ) {
          expl[nm] = exs[nm] + exn[nm];
          exmi[nm] = exs[nm] - exn[nm];
        }
      }
      ntot = 0;
      for( nl=0; nl<mp; nl++ ) {
        for( nm=0; nm<mp; nm+=2 ) {
          mmax = nfour[nl]-1;
          if( mmax > nm ) mmax = nm;
          for( mth=0; mth<=mmax; mth++ ) {
            nms = nm*(nm+1)/2+mth;
            ncurrent = ntot+mth;
            aaxd[nms] += m2e[nl][mth][nm]*wt[nl]*expl[ncurrent];
          }
        }
        for( nm=1; nm<mp; nm+=2 ) {
          mmax = nfour[nl]-1;
          if( mmax > nm ) mmax = nm;
          for( mth=0; mth<=mmax; mth++ ) {
            nms = nm*(nm+1)/2+mth;
            ncurrent = ntot+mth;
            aaxd[nms] += m2e[nl][mth][nm]*wt[nl]*exmi[ncurrent];
          }
        }
        ntot += nfour[nl];
      }
      rscale = pi/lscale*0.5;
      for( nm=0; nm<mp; nm++ ) {
        for( mth=0; mth<=nm; mth++ ) {
          nms = nm*(nm+1)/2+mth;
          aaxd[nms] *= zeye[mth]*rscale*double(ytop[nm][mth]);
        }
      }
    }

// ROTYTOZ
    for( l=0; l<mp; l++ ) {
      for( m=0; m<=l; m++ ) {
        ln = l*(l+1)/2;
        lm = l*(l+1)/2+m;
        aax[lm] = aaxd[ln]*rdpi2[l][0][mp+m];
        for( n=1; n<=l; n++ ) {
          ln = l*(l+1)/2+n;
          aax[lm] += aaxd[ln]*rdpi2[l][n][mp+m]+std::conj(aaxd[ln])*rdpi2[l][n][mp-m];
        }
      }
    }
    ephi[0] = 1;
    for( m=1; m<mp+1; m++ ) {
      ephi[m] = ephi[m-1]*ima;
    }
    for( l=0; l<mp; l++ ) {
      for( m=0; m<=l; m++ ) {
        lm = l*(l+1)/2+m;
        aaxd[lm] = ephi[m]*aax[lm];
        ax[ii][lm] += aaxd[lm];
      }
    }
  }

// process east and west lists ---------------------------------------------------------------------
  for( ii=0; ii<lbi; ii++ ) {
    ib = ii+nlbi[lev-1];
    iexp1[ib] = 0;
    iexp2[ib] = 0;
    for( jj=0; jj<maxwave; jj++ ) {
      lexp1[ib][jj] = 0;
      lexp2[ib][jj] = 0;
    }
  }

  nbox = -1;
  for( ii=0; ii<lbi; ii++ ) {
    ib = ii+nlbi[lev-1];
    for( i=0; i<nmp; i++ ) {
      bbxd[i] = bx[ib][i];
    }

// ROTZTOX
    for( l=0; l<mp; l++ ) {
      for( m=0; m<=l; m++ ) {
        ln = l*(l+1)/2;
        lm = l*(l+1)/2+m;
        bbx[lm] = bbxd[ln]*rdpi2[l][0][mp+m];
        for( n=1; n<=l; n++ ) {
          ln = l*(l+1)/2+n;
          bbx[lm] += bbxd[ln]*rdpi2[l][n][mp+m]+std::conj(bbxd[ln])*rdpi2[l][n][mp-m];
        }
      }
    }

// YMKEWEXP
// YMPOLETOEXP
    ntot = 0;
    for( l=0; l<mp; l++ ) {
      sgn = -1;
      for( m=0; m<nfour[l]; m++ ) {
        ic = ntot+m;
        ztmp1 = 0;
        ztmp2 = 0;
        sgn = -sgn;
        for( n=m; n<mp; n+=2 ) {
          nm = n*(n+1)/2+m;
          ztmp1 += std::complex<double>(m2e[l][m][n])*bbx[nm];
        }
        for( n=m+1; n<mp; n+=2 ) {
          nm = n*(n+1)/2+m;
          ztmp2 += std::complex<double>(m2e[l][m][n])*bbx[nm];
        }
        exe[ic] = ztmp1+ztmp2;
        exw[ic] = sgn*(ztmp1-ztmp2);
      }
      ntot += nfour[l];
    }

// YFTOPHYS
    nftot = 0;
    nptot = 0;
    nexte = 0;
    nexto = 0;
    for( l=0; l<mp; l++ ) {
      for( ival=0; ival<nwave[l]/2; ival++ ) {
        fxe[nptot+ival] = exe[nftot];
        fxw[nptot+ival] = exw[nftot];
        sgn = -2;
        for( mm=1; mm<nfour[l]; mm+=2 ) {
          sgn = -sgn;
          fxe[nptot+ival] += ima*std::complex<double>(sgn*std::real( fexpe[nexte]*exe[nftot+mm] ));
          fxw[nptot+ival] += ima*std::complex<double>(sgn*std::real( fexpe[nexte]*exw[nftot+mm] ));
          nexte++;
        }
        sgn = 2;
        for( mm=2; mm<nfour[l]; mm+=2 ) {
          sgn = -sgn;
          fxe[nptot+ival] += sgn*std::real( fexpo[nexto]*exe[nftot+mm] );;
          fxw[nptot+ival] += sgn*std::real( fexpo[nexto]*exw[nftot+mm] );;
          nexto++;
        }
      }
      nftot += nfour[l];
      nptot += nwave[l]/2;
    }

// YMKEWEXP part2
    islast = 1;
    if( nbox != nfi[ib]/8 ) {
      nbox = nfi[ib]/8;
      iexp[0]=0;
      iexp[1]=0;
      iexp[2]=0;
      iexp[3]=0;
      iexp[4]=0;
      iexp[5]=0;
      iexp[6]=0;
      iexp[7]=0;
      iexp[8]=0;
      iexp[9]=0;
      iexp[10]=0;
      iexp[11]=0;
      iexp[12]=0;
      iexp[13]=0;
      iexp[14]=0;
      iexp[15]=0;
      for( jj=0; jj<maxwave; jj++ ) {
        fxeall[jj] = 0;
        fxe1357[jj] = 0;
        fxe13[jj] = 0;
        fxe57[jj] = 0;
        fxe1[jj] = 0;
        fxe3[jj] = 0;
        fxe5[jj] = 0;
        fxe7[jj] = 0;
        fxwall[jj] = 0;
        fxw2468[jj] = 0;
        fxw24[jj] = 0;
        fxw68[jj] = 0;
        fxw2[jj] = 0;
        fxw4[jj] = 0;
        fxw6[jj] = 0;
        fxw8[jj] = 0;
      }
    }

    if( nfi[ib] % 8 == 0 ) {
      for( jj=0; jj<maxwave; jj++ ) {
        fxe1[jj] = fxe[jj];
        fxwall[jj] += fxw[jj];
      }
      iexp[4]++;
      iexp[8]++;
    } else if( nfi[ib] % 8 == 1 ) {
      for( jj=0; jj<maxwave; jj++ ) {
        fxeall[jj] += fxe[jj]/zs[jj][0];
        fxw2[jj] = fxw[jj]*zs[jj][0];
      }
      iexp[0]++;
      iexp[12]++;
    } else if( nfi[ib] % 8 == 2 ) {
      for( jj=0; jj<maxwave; jj++ ) {
        fxe3[jj] = fxe[jj]*std::conj(ys[jj][0]);
        fxwall[jj] += fxw[jj]*ys[jj][0];
      }
      iexp[5]++;
      iexp[8]++;
    } else if( nfi[ib] % 8 == 3 ) {
      for( jj=0; jj<maxwave; jj++ ) {
        ztmp = ys[jj][0]*zs[jj][0];
        fxeall[jj] += fxe[jj]/ztmp;
        fxw4[jj] = fxw[jj]*ztmp;
      }
      iexp[0]++;
      iexp[13]++;
    } else if( nfi[ib] % 8 == 4 ) {
      for( jj=0; jj<maxwave; jj++ ) {
        fxe5[jj] = fxe[jj]*xs[jj][0];
        fxwall[jj] += fxw[jj]*std::conj(xs[jj][0]);
      }
      iexp[6]++;
      iexp[8]++;
    } else if( nfi[ib] % 8 == 5 ) {
      for( jj=0; jj<maxwave; jj++ ) {
        ztmp = xs[jj][0]/zs[jj][0];
        fxeall[jj] += fxe[jj]*ztmp;
        fxw6[jj] = fxw[jj]/ztmp;
      }
      iexp[0]++;
      iexp[14]++;
    } else if( nfi[ib] % 8 == 6 ) {
      for( jj=0; jj<maxwave; jj++ ) {
        ztmp = xs[jj][0]*std::conj(ys[jj][0]);
        fxe7[jj] = fxe[jj]*ztmp;
        fxwall[jj] += fxw[jj]*std::conj(ztmp);
      }
      iexp[7]++;
      iexp[8]++;
    } else if( nfi[ib] % 8 == 7 ) {
      for( jj=0; jj<maxwave; jj++ ) {
        ztmp = xs[jj][0]*std::conj(ys[jj][0])/zs[jj][0];
        fxeall[jj] += fxe[jj]*ztmp;
        fxw8[jj] = fxw[jj]/ztmp;
      }
      iexp[15]++;
    }
    if( ib == lbi-1 ) {
      for( jj=0; jj<maxwave; jj++ ) {
        fxe13[jj] = fxe1[jj]+fxe3[jj];
        fxe57[jj] = fxe5[jj]+fxe7[jj];
        fxe1357[jj] = fxe13[jj]+fxe57[jj];
        fxeall[jj] += fxe1357[jj];
        fxw24[jj] = fxw2[jj]+fxw4[jj];
        fxw68[jj] = fxw6[jj]+fxw8[jj];
        fxw2468[jj] = fxw24[jj]+fxw68[jj];
        fxwall[jj] += fxw2468[jj];
      }
      islast = 0;
      iexp[0]++;
      iexp[1]++;
      iexp[2]++;
      iexp[3]++;
      iexp[8]++;
      iexp[9]++;
      iexp[10]++;
      iexp[11]++;
    } else if ( nfi[ib+1]/8 != nbox ) {
      for( jj=0; jj<maxwave; jj++ ) {
        fxe13[jj] = fxe1[jj]+fxe3[jj];
        fxe57[jj] = fxe5[jj]+fxe7[jj];
        fxe1357[jj] = fxe13[jj]+fxe57[jj];
        fxeall[jj] += fxe1357[jj];
        fxw24[jj] = fxw2[jj]+fxw4[jj];
        fxw68[jj] = fxw6[jj]+fxw8[jj];
        fxw2468[jj] = fxw24[jj]+fxw68[jj];
        fxwall[jj] += fxw2468[jj];
      }
      islast = 0;
      iexp[0]++;
      iexp[1]++;
      iexp[2]++;
      iexp[3]++;
      iexp[8]++;
      iexp[9]++;
      iexp[10]++;
      iexp[11]++;
    }

    if( islast == 0 ) {

// MKEALIST
      neall = 0;
      ne1357 = 0;
      ne13 = 0;
      ne57 = 0;
      ne1 = 0;
      ne3 = 0;
      ne5 = 0;
      ne7 = 0;
      boxc(nfi[ib]/8,3,nc);
      ixp = nc[0]+1;
      iyp = nc[1]+1;
      izp = nc[2]+1;
      jxp = ixp+1;
      for( jzp=izp-1; jzp<=izp+1; jzp++ ) {
        for( jyp=iyp-1; jyp<=iyp+1; jyp++ ) {
          for( jz=std::max(2*jzp-2,jzmin); jz<=std::min(2*jzp-1,jzmax); jz++ ) {
            for( jy=std::max(2*jyp-2,jymin); jy<=std::min(2*jyp-1,jymax); jy++ ) {
              for( jx=std::max(2*jxp-2,jxmin); jx<=std::min(2*jxp-1,jxmax); jx++ ) {
                ix = jx-2*(ixp-1);
                iy = jy-2*(iyp-1);
                iz = jz-2*(izp-1);
                if( ix > 2.5 ) {
                  if( iz >= -1 && iz <= 2 && iy >= -1 && iy <= 2 ) {
                    nc[0] = jx;
                    nc[1] = jy;
                    nc[2] = jz;
                    boxn1(nc,je,lev);
                    jj = nej[je];
                    if( jj != -1 ) {
                      ieall[neall] = jj;
                      ixeall[neall] = -iz;
                      iyeall[neall] = iy;
                      neall++;
                    }
                  }
                } else if( ( iz == 0 || iz == 1 ) && ( iy == 0 || iy == 1 ) ) {
                  nc[0] = jx;
                  nc[1] = jy;
                  nc[2] = jz;
                  boxn1(nc,je,lev);
                  jj = nej[je];
                  if( jj != -1 ) {
                    ie1357[ne1357] = jj;
                    ixe1357[ne1357] = -iz;
                    iye1357[ne1357] = iy;
                    ne1357++;
                  }
                } else if( iz == -1 && ( iy == 0 || iy == 1 ) ) {
                  nc[0] = jx;
                  nc[1] = jy;
                  nc[2] = jz;
                  boxn1(nc,je,lev);
                  jj = nej[je];
                  if( jj != -1 ) {
                    ie13[ne13] = jj;
                    ixe13[ne13] = -iz;
                    iye13[ne13] = iy;
                    ne13++;
                  }
                } else if( iz == 2 && ( iy == 0 || iy == 1 ) ) {
                  nc[0] = jx;
                  nc[1] = jy;
                  nc[2] = jz;
                  boxn1(nc,je,lev);
                  jj = nej[je];
                  if( jj != -1 ) {
                    ie57[ne57] = jj;
                    ixe57[ne57] = -iz;
                    iye57[ne57] = iy;
                    ne57++;
                  }
                } else if( iy == -1 ) {
                  if( iz >= -1 && iz <= 1 ) {
                    nc[0] = jx;
                    nc[1] = jy;
                    nc[2] = jz;
                    boxn1(nc,je,lev);
                    jj = nej[je];
                    if( jj != -1 ) {
                      ie1[ne1] = jj;
                      ixe1[ne1] = -iz;
                      iye1[ne1] = iy;
                      ne1++;
                    }
                  }
                  if( iz >= 0 && iz <= 2 ) {
                    nc[0] = jx;
                    nc[1] = jy;
                    nc[2] = jz;
                    boxn1(nc,je,lev);
                    jj = nej[je];
                    if( jj != -1 ) {
                      ie5[ne5] = jj;
                      ixe5[ne5] = -iz;
                      iye5[ne5] = iy;
                      ne5++;
                    }
                  }
                } else if( iy == 2 ) {
                  if( iz >= -1 && iz <= 1 ) {
                    nc[0] = jx;
                    nc[1] = jy;
                    nc[2] = jz;
                    boxn1(nc,je,lev);
                    jj = nej[je];
                    if( jj != -1 ) {
                      ie3[ne3] = jj;
                      ixe3[ne3] = -iz;
                      iye3[ne3] = iy;
                      ne3++;
                    }
                  }
                  if( iz >= 0 && iz <= 2 ) {
                    nc[0] = jx;
                    nc[1] = jy;
                    nc[2] = jz;
                    boxn1(nc,je,lev);
                    jj = nej[je];
                    if( jj != -1 ) {
                      ie7[ne7] = jj;
                      ixe7[ne7] = -iz;
                      iye7[ne7] = iy;
                      ne7++;
                    }
                  }
                }
              }
            }
          }
        }
      }

// YPROCESSEA
      if( iexp[0] > 0 ) {
        for( i=0; i<neall; i++ ) {
          iexp1[ieall[i]]++;
          for( jj=0; jj<maxwave; jj++ ) {
            zmul = zs[jj][2];
            if( ixeall[i] > 0 ) zmul *= xs[jj][ixeall[i]-1];
            if( ixeall[i] < 0 ) zmul *= std::conj(xs[jj][-ixeall[i]-1]);
            if( iyeall[i] > 0 ) zmul *= ys[jj][iyeall[i]-1];
            if( iyeall[i] < 0 ) zmul *= std::conj(ys[jj][-iyeall[i]-1]);
            lexp1[ieall[i]][jj] += fxeall[jj]*zmul;
          }
        }
      }
      if( iexp[1] > 0 ) {
        for( i=0; i<ne1357; i++ ) {
          iexp1[ie1357[i]]++;
          for( jj=0; jj<maxwave; jj++ ) {
            zmul = zs[jj][1];
            if( ixe1357[i] > 0 ) zmul *= xs[jj][ixe1357[i]-1];
            if( ixe1357[i] < 0 ) zmul *= std::conj(xs[jj][-ixe1357[i]-1]);
            if( iye1357[i] > 0 ) zmul *= ys[jj][iye1357[i]-1];
            if( iye1357[i] < 0 ) zmul *= std::conj(ys[jj][-iye1357[i]-1]);
            lexp1[ie1357[i]][jj] += fxe1357[jj]*zmul;
          }
        }
      }
      if( iexp[2] > 0 ) {
        for( i=0; i<ne13; i++ ) {
          iexp1[ie13[i]]++;
          for( jj=0; jj<maxwave; jj++ ) {
            zmul = zs[jj][1];
            if( ixe13[i] > 0 ) zmul *= xs[jj][ixe13[i]-1];
            if( ixe13[i] < 0 ) zmul *= std::conj(xs[jj][-ixe13[i]-1]);
            if( iye13[i] > 0 ) zmul *= ys[jj][iye13[i]-1];
            if( iye13[i] < 0 ) zmul *= std::conj(ys[jj][-iye13[i]-1]);
            lexp1[ie13[i]][jj] += fxe13[jj]*zmul;
          }
        }
      }
      if( iexp[3] > 0 ) {
        for( i=0; i<ne57; i++ ) {
          iexp1[ie57[i]]++;
          for( jj=0; jj<maxwave; jj++ ) {
            zmul = zs[jj][1];
            if( ixe57[i] > 0 ) zmul *= xs[jj][ixe57[i]-1];
            if( ixe57[i] < 0 ) zmul *= std::conj(xs[jj][-ixe57[i]-1]);
            if( iye57[i] > 0 ) zmul *= ys[jj][iye57[i]-1];
            if( iye57[i] < 0 ) zmul *= std::conj(ys[jj][-iye57[i]-1]);
            lexp1[ie57[i]][jj] += fxe57[jj]*zmul;
          }
        }
      }
      if( iexp[4] > 0 ) {
        for( i=0; i<ne1; i++ ) {
          iexp1[ie1[i]]++;
          for( jj=0; jj<maxwave; jj++ ) {
            zmul = zs[jj][1];
            if( ixe1[i] > 0 ) zmul *= xs[jj][ixe1[i]-1];
            if( ixe1[i] < 0 ) zmul *= std::conj(xs[jj][-ixe1[i]-1]);
            if( iye1[i] > 0 ) zmul *= ys[jj][iye1[i]-1];
            if( iye1[i] < 0 ) zmul *= std::conj(ys[jj][-iye1[i]-1]);
            lexp1[ie1[i]][jj] += fxe1[jj]*zmul;
          }
        }
      }
      if( iexp[5] > 0 ) {
        for( i=0; i<ne3; i++ ) {
          iexp1[ie3[i]]++;
          for( jj=0; jj<maxwave; jj++ ) {
            zmul = zs[jj][1];
            if( ixe3[i] > 0 ) zmul *= xs[jj][ixe3[i]-1];
            if( ixe3[i] < 0 ) zmul *= std::conj(xs[jj][-ixe3[i]-1]);
            if( iye3[i] > 0 ) zmul *= ys[jj][iye3[i]-1];
            if( iye3[i] < 0 ) zmul *= std::conj(ys[jj][-iye3[i]-1]);
            lexp1[ie3[i]][jj] += fxe3[jj]*zmul;
          }
        }
      }
      if( iexp[6] > 0 ) {
        for( i=0; i<ne5; i++ ) {
          iexp1[ie5[i]]++;
          for( jj=0; jj<maxwave; jj++ ) {
            zmul = zs[jj][1];
            if( ixe5[i] > 0 ) zmul *= xs[jj][ixe5[i]-1];
            if( ixe5[i] < 0 ) zmul *= std::conj(xs[jj][-ixe5[i]-1]);
            if( iye5[i] > 0 ) zmul *= ys[jj][iye5[i]-1];
            if( iye5[i] < 0 ) zmul *= std::conj(ys[jj][-iye5[i]-1]);
            lexp1[ie5[i]][jj] += fxe5[jj]*zmul;
          }
        }
      }
      if( iexp[7] > 0 ) {
        for( i=0; i<ne7; i++ ) {
          iexp1[ie7[i]]++;
          for( jj=0; jj<maxwave; jj++ ) {
            zmul = zs[jj][1];
            if( ixe7[i] > 0 ) zmul *= xs[jj][ixe7[i]-1];
            if( ixe7[i] < 0 ) zmul *= std::conj(xs[jj][-ixe7[i]-1]);
            if( iye7[i] > 0 ) zmul *= ys[jj][iye7[i]-1];
            if( iye7[i] < 0 ) zmul *= std::conj(ys[jj][-iye7[i]-1]);
            lexp1[ie7[i]][jj] += fxe7[jj]*zmul;
          }
        }
      }

// MKWELIST
      nwall = 0;
      nw2468 = 0;
      nw24 = 0;
      nw68 = 0;
      nw2 = 0;
      nw4 = 0;
      nw6 = 0;
      nw8 = 0;
      boxc(nfi[ib]/8,3,nc);
      ixp = nc[0]+1;
      iyp = nc[1]+1;
      izp = nc[2]+1;
      jxp = ixp-1;
      for( jzp=izp-1; jzp<=izp+1; jzp++ ) {
        for( jyp=iyp-1; jyp<=iyp+1; jyp++ ) {
          for( jz=std::max(2*jzp-2,jzmin); jz<=std::min(2*jzp-1,jzmax); jz++ ) {
            for( jy=std::max(2*jyp-2,jymin); jy<=std::min(2*jyp-1,jymax); jy++ ) {
              for( jx=std::max(2*jxp-2,jxmin); jx<=std::min(2*jxp-1,jxmax); jx++ ) {
                ix = jx-2*(ixp-1);
                iy = jy-2*(iyp-1);
                iz = jz-2*(izp-1);
                if( ix <= -1.5 ) {
                  if( iz >= -1 && iz <= 2 && iy >= -1 && iy <= 2 ) {
                    nc[0] = jx;
                    nc[1] = jy;
                    nc[2] = jz;
                    boxn1(nc,je,lev);
                    jj = nej[je];
                    if( jj != -1 ) {
                      iwall[nwall] = jj;
                      ixwall[nwall] = -iz;
                      iywall[nwall] = iy;
                      nwall++;
                    }
                  }
                } else if( ( iz == 0 || iz == 1 ) && ( iy == 0 || iy == 1 ) ){
                  nc[0] = jx;
                  nc[1] = jy;
                  nc[2] = jz;
                  boxn1(nc,je,lev);
                  jj = nej[je];
                  if( jj != -1 ) {
                    iw2468[nw2468] = jj;
                    ixw2468[nw2468] = -iz;
                    iyw2468[nw2468] = iy;
                    nw2468++;
                  }
                } else if( iz == -1 && ( iy == 0 || iy == 1 ) ) {
                  nc[0] = jx;
                  nc[1] = jy;
                  nc[2] = jz;
                  boxn1(nc,je,lev);
                  jj = nej[je];
                  if( jj != -1 ) {
                    iw24[nw24] = jj;
                    ixw24[nw24] = -iz;
                    iyw24[nw24] = iy;
                    nw24++;
                  }
                } else if( iz == 2 && ( iy == 0 || iy == 1 ) ) {
                  nc[0] = jx;
                  nc[1] = jy;
                  nc[2] = jz;
                  boxn1(nc,je,lev);
                  jj = nej[je];
                  if( jj != -1 ) {
                    iw68[nw68] = jj;
                    ixw68[nw68] = -iz;
                    iyw68[nw68] = iy;
                    nw68++;
                  }
                } else if( iy == -1 ) {
                  if( iz >= -1 && iz <= 1 ) {
                    nc[0] = jx;
                    nc[1] = jy;
                    nc[2] = jz;
                    boxn1(nc,je,lev);
                    jj = nej[je];
                    if( jj != -1 ) {
                      iw2[nw2] = jj;
                      ixw2[nw2] = -iz;
                      iyw2[nw2] = iy;
                      nw2++;
                    }
                  }
                  if( iz >= 0 && iz <= 2 ) {
                    nc[0] = jx;
                    nc[1] = jy;
                    nc[2] = jz;
                    boxn1(nc,je,lev);
                    jj = nej[je];
                    if( jj != -1 ) {
                      iw6[nw6] = jj;
                      ixw6[nw6] = -iz;
                      iyw6[nw6] = iy;
                      nw6++;
                    }
                  }
                } else if( iy == 2 ) {
                  if( iz >= -1 && iz <= 1 ) {
                    nc[0] = jx;
                    nc[1] = jy;
                    nc[2] = jz;
                    boxn1(nc,je,lev);
                    jj = nej[je];
                    if( jj != -1 ) {
                      iw4[nw4] = jj;
                      ixw4[nw4] = -iz;
                      iyw4[nw4] = iy;
                      nw4++;
                    }
                  }
                  if( iz >= 0 && iz <= 2 ) {
                    nc[0] = jx;
                    nc[1] = jy;
                    nc[2] = jz;
                    boxn1(nc,je,lev);
                    jj = nej[je];
                    if( jj != -1 ) {
                      iw8[nw8] = jj;
                      ixw8[nw8] = -iz;
                      iyw8[nw8] = iy;
                      nw8++;
                    }
                  }
                }
              }
            }
          }
        }
      }

// YPROCESSWE
      if( iexp[8] > 0 ) {
        for( i=0; i<nwall; i++ ) {
          iexp2[iwall[i]]++;
          for( jj=0; jj<maxwave; jj++ ) {
            zmul = zs[jj][1];
            if( ixwall[i] > 0 ) zmul *= std::conj(xs[jj][ixwall[i]-1]);
            if( ixwall[i] < 0 ) zmul *= xs[jj][-ixwall[i]-1];
            if( iywall[i] > 0 ) zmul *= std::conj(ys[jj][iywall[i]-1]);
            if( iywall[i] < 0 ) zmul *= ys[jj][-iywall[i]-1];
            lexp2[iwall[i]][jj] += fxwall[jj]*zmul;
          }
        }
      }
      if( iexp[9] > 0 ) {
        for( i=0; i<nw2468; i++ ) {
          iexp2[iw2468[i]]++;
          for( jj=0; jj<maxwave; jj++ ) {
            zmul = zs[jj][0];
            if( ixw2468[i] > 0 ) zmul *= std::conj(xs[jj][ixw2468[i]-1]);
            if( ixw2468[i] < 0 ) zmul *= xs[jj][-ixw2468[i]-1];
            if( iyw2468[i] > 0 ) zmul *= std::conj(ys[jj][iyw2468[i]-1]);
            if( iyw2468[i] < 0 ) zmul *= ys[jj][-iyw2468[i]-1];
            lexp2[iw2468[i]][jj] += fxw2468[jj]*zmul;
          }
        }
      }
      if( iexp[10] > 0 ) {
        for( i=0; i<nw24; i++ ) {
          iexp2[iw24[i]]++;
          for( jj=0; jj<maxwave; jj++ ) {
            zmul = zs[jj][0];
            if( ixw24[i] > 0 ) zmul *= std::conj(xs[jj][ixw24[i]-1]);
            if( ixw24[i] < 0 ) zmul *= xs[jj][-ixw24[i]-1];
            if( iyw24[i] > 0 ) zmul *= std::conj(ys[jj][iyw24[i]-1]);
            if( iyw24[i] < 0 ) zmul *= ys[jj][-iyw24[i]-1];
            lexp2[iw24[i]][jj] += fxw24[jj]*zmul;
          }
        }
      }
      if( iexp[11] > 0 ) {
        for( i=0; i<nw68; i++ ) {
          iexp2[iw68[i]]++;
          for( jj=0; jj<maxwave; jj++ ) {
            zmul = zs[jj][0];
            if( ixw68[i] > 0 ) zmul *= std::conj(xs[jj][ixw68[i]-1]);
            if( ixw68[i] < 0 ) zmul *= xs[jj][-ixw68[i]-1];
            if( iyw68[i] > 0 ) zmul *= std::conj(ys[jj][iyw68[i]-1]);
            if( iyw68[i] < 0 ) zmul *= ys[jj][-iyw68[i]-1];
            lexp2[iw68[i]][jj] += fxw68[jj]*zmul;
          }
        }
      }
      if( iexp[12] > 0 ) {
        for( i=0; i<nw2; i++ ) {
          iexp2[iw2[i]]++;
          for( jj=0; jj<maxwave; jj++ ) {
            zmul = zs[jj][0];
            if( ixw2[i] > 0 ) zmul *= std::conj(xs[jj][ixw2[i]-1]);
            if( ixw2[i] < 0 ) zmul *= xs[jj][-ixw2[i]-1];
            if( iyw2[i] > 0 ) zmul *= std::conj(ys[jj][iyw2[i]-1]);
            if( iyw2[i] < 0 ) zmul *= ys[jj][-iyw2[i]-1];
            lexp2[iw2[i]][jj] += fxw2[jj]*zmul;
          }
        }
      }
      if( iexp[13] > 0 ) {
        for( i=0; i<nw4; i++ ) {
          iexp2[iw4[i]]++;
          for( jj=0; jj<maxwave; jj++ ) {
            zmul = zs[jj][0];
            if( ixw4[i] > 0 ) zmul *= std::conj(xs[jj][ixw4[i]-1]);
            if( ixw4[i] < 0 ) zmul *= xs[jj][-ixw4[i]-1];
            if( iyw4[i] > 0 ) zmul *= std::conj(ys[jj][iyw4[i]-1]);
            if( iyw4[i] < 0 ) zmul *= ys[jj][-iyw4[i]-1];
            lexp2[iw4[i]][jj] += fxw4[jj]*zmul;
          }
        }
      }
      if( iexp[14] > 0 ) {
        for( i=0; i<nw6; i++ ) {
          iexp2[iw6[i]]++;
          for( jj=0; jj<maxwave; jj++ ) {
            zmul = zs[jj][0];
            if( ixw6[i] > 0 ) zmul *= std::conj(xs[jj][ixw6[i]-1]);
            if( ixw6[i] < 0 ) zmul *= xs[jj][-ixw6[i]-1];
            if( iyw6[i] > 0 ) zmul *= std::conj(ys[jj][iyw6[i]-1]);
            if( iyw6[i] < 0 ) zmul *= ys[jj][-iyw6[i]-1];
            lexp2[iw6[i]][jj] += fxw6[jj]*zmul;
          }
        }
      }
      if( iexp[15] > 0 ) {
        for( i=0; i<nw8; i++ ) {
          iexp2[iw8[i]]++;
          for( jj=0; jj<maxwave; jj++ ) {
            zmul = zs[jj][0];
            if( ixw8[i] > 0 ) zmul *= std::conj(xs[jj][ixw8[i]-1]);
            if( ixw8[i] < 0 ) zmul *= xs[jj][-ixw8[i]-1];
            if( iyw8[i] > 0 ) zmul *= std::conj(ys[jj][iyw8[i]-1]);
            if( iyw8[i] < 0 ) zmul *= ys[jj][-iyw8[i]-1];
            lexp2[iw8[i]][jj] += fxw8[jj]*zmul;
          }
        }
      }
    }
  }

  for( ii=0; ii<lbi; ii++ ) {
    if( iexp1[ii] > 0 ) {
// YPHYSTOF
      nftot = 0;
      nptot = 0;
      next = 0;
      for( l=0; l<mp; l++ ) {
        nalpha = nwave[l];
        nalpha2 = nalpha/2;
        halpha = 2*pi/nalpha;
        exe[nftot] = 0;
        for( ival=0; ival<nalpha2; ival++ ) {
          exe[nftot] += 2.0*std::real(lexp1[ii][nptot+ival]);
        }
        exe[nftot] /= nalpha;
        for( mm=2; mm<nfour[l]; mm+=2 ) {
          exe[nftot+mm] = 0;
          for( ival=0; ival<nalpha2; ival++ ) {
            exe[nftot+mm] += fexpb[next]*2.0*std::real(lexp1[ii][nptot+ival]);
            next++;
          }
          exe[nftot+mm] /= nalpha;
        }
        for( mm=1; mm<nfour[l]; mm+=2 ) {
          exe[nftot+mm] = 0;
          for( ival=0; ival<nalpha2; ival++ ) {
            exe[nftot+mm] += fexpb[next]*2.0*ima*std::complex<double>(std::imag(lexp1[ii][nptot+ival]));
            next++;
          }
          exe[nftot+mm] /= nalpha;
        }
        nftot += nfour[l];
        nptot += nwave[l]/2;
      }
    }
    if( iexp2[ii] > 0 ) {
// YPHYSTOF
      nftot = 0;
      nptot = 0;
      next = 0;
      for( l=0; l<mp; l++ ) {
        nalpha = nwave[l];
        nalpha2 = nalpha/2;
        halpha = 2*pi/nalpha;
        exw[nftot] = 0;
        for( ival=0; ival<nalpha2; ival++ ) {
          exw[nftot] += 2.0*std::real(lexp2[ii][nptot+ival]);
        }
        exw[nftot] /= nalpha;
        for( mm=2; mm<nfour[l]; mm+=2 ) {
          exw[nftot+mm] = 0;
          for( ival=0; ival<nalpha2; ival++ ) {
            exw[nftot+mm] += fexpb[next]*2.0*std::real(lexp2[ii][nptot+ival]);
            next++;
          }
          exw[nftot+mm] /= nalpha;
        }
        for( mm=1; mm<nfour[l]; mm+=2 ) {
          exw[nftot+mm] = 0;
          for( ival=0; ival<nalpha2; ival++ ) {
            exw[nftot+mm] += fexpb[next]*2.0*ima*std::complex<double>(std::imag(lexp2[ii][nptot+ival]));
            next++;
          }
          exw[nftot+mm] /= nalpha;
        }
        nftot += nfour[l];
        nptot += nwave[l]/2;
      }
    }
    if( iexp1[ii] > 0 || iexp2[ii] > 0 ) {
// YEXPTOLOCAL
      zeye[0] = 1;
      for( i=1; i<nmax; i++ ) {
        zeye[i] = zeye[i-1]*ima;
      }
      for( nm=0; nm<mp; nm++ ) {
        for( mth=0; mth<mp; mth++ ) {
          nms = nm*(nm+1)/2+mth;
          aaxd[nms] = 0;
        }
      }
      if( iexp1[ii] <= 0 ) {
        for( nm=0; nm<maxfour; nm++ ) {
          expl[nm] = exw[nm];
          exmi[nm] = exw[nm];
        }
      } else if( iexp2[ii] <= 0 ) {
        for( nm=0; nm<maxfour; nm++ ) {
          expl[nm] = exe[nm];
          exmi[nm] = -exe[nm];
        }
      } else {
        for( nm=0; nm<maxfour; nm++ ) {
          expl[nm] = exw[nm] + exe[nm];
          exmi[nm] = exw[nm] - exe[nm];
        }
      }
      ntot = 0;
      for( nl=0; nl<mp; nl++ ) {
        for( nm=0; nm<mp; nm+=2 ) {
          mmax = nfour[nl]-1;
          if( mmax > nm ) mmax = nm;
          for( mth=0; mth<=mmax; mth++ ) {
            nms = nm*(nm+1)/2+mth;
            ncurrent = ntot+mth;
            aaxd[nms] += m2e[nl][mth][nm]*wt[nl]*expl[ncurrent];
          }
        }
        for( nm=1; nm<mp; nm+=2 ) {
          mmax = nfour[nl]-1;
          if( mmax > nm ) mmax = nm;
          for( mth=0; mth<=mmax; mth++ ) {
            nms = nm*(nm+1)/2+mth;
            ncurrent = ntot+mth;
            aaxd[nms] += m2e[nl][mth][nm]*wt[nl]*exmi[ncurrent];
          }
        }
        ntot += nfour[nl];
      }
      rscale = pi/lscale*0.5;
      for( nm=0; nm<mp; nm++ ) {
        for( mth=0; mth<=nm; mth++ ) {
          nms = nm*(nm+1)/2+mth;
          aaxd[nms] *= zeye[mth]*rscale*double(ytop[nm][mth]);
        }
      }
    }

// ROTYTOZ
    for( l=0; l<mp; l++ ) {
      for( m=0; m<=l; m++ ) {
        ln = l*(l+1)/2;
        lm = l*(l+1)/2+m;
        aax[lm] = aaxd[ln]*rdmpi2[l][0][mp+m];
        for( n=1; n<=l; n++ ) {
          ln = l*(l+1)/2+n;
          aax[lm] += aaxd[ln]*rdmpi2[l][n][mp+m]+std::conj(aaxd[ln])*rdmpi2[l][n][mp-m];
        }
        ax[ii][lm] += aax[lm];
      }
    }
  }

  for( jj=0; jj<lbj; jj++ ) {
    jb = njb[jj];
    for( j=0; j<nmp; j++ ) {
      bx[jb][j] = 0;
    }
  }
  delete[] fexpe;
  delete[] fexpo;
  delete[] fexpb;
  delete[] xs;
  delete[] ys;
  delete[] zs;
  for( i=0; i<nbne; i++ ) delete[] lexp1[i];
  delete[] lexp1;
  for( i=0; i<nbne; i++ ) delete[] lexp2[i];
  delete[] lexp2;

  toc = get_time();
  tfmm[20] += toc-tic;
}
