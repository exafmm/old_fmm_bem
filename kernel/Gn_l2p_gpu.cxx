#include "../include/constants.h"

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi,*vi;
extern int **ndi,*nfi,*nc;
extern std::complex<double> (*ax)[mpsym];
extern float *fac;
extern double *tfmm;
extern int *istagp,*iendgp,*nvecd;
extern float *xig,*yig,*zig,*gxig,*gyig,*gzig,*vig;
extern float *arex,*aimx;

extern void l2pgpu_(int*, double*, double*, double*, double*, double*,
        float*, float*, float*, float*, float*, float*, float*,
        float*, float*, float*);

void Gn_l2p(int nmp, int mp, int lbi) {
  int kl,idev,mblok,ni,nj,ncall,ii,icall,iblok,jc,jbase,j,jsize,ibase,isize,is,i;
  double rb,op,xmind,ymind,zmind;
  tic = get_time();

  kl = 1 << lmax;
  rb = rd/kl;
  idev = myrank%ngpu;
  mblok = 3;

  ni = 0;
  nj = 0;
  ncall = 0;
  istagp[0] = 0;
  for( ii=0; ii<lbi; ii++ ) {
    ni += ((ndi[1][ii]-ndi[0][ii]+nblok1)/nblok1+1)*nblok1;
    nj += nmp;
    if( ni > nimax || nj > njmax ) {
      iendgp[ncall] = ii-1;
      ncall++;
      istagp[ncall] = ii;
      ni = ((ndi[1][ii]-ndi[0][ii]+nblok1)/nblok1+1)*nblok1;
      nj = nmp;
    }
  }
  iendgp[ncall] = lbi-1;
  if(lbi != 0) ncall++;

  for( icall=0; icall<ncall; icall++ ) {
    iblok = 0;
    jc = 0;
    op = 0;
    for( ii=istagp[icall]; ii<=iendgp[icall]; ii++ ) {
      jbase = jc;
      for( j=0; j<nmp; j++ ) {
        arex[jc] = std::real(ax[ii][j]);
        aimx[jc] = std::imag(ax[ii][j]);
        jc++;
      }
      jsize = jc-jbase;
      ibase = ndi[0][ii];
      isize = ndi[1][ii]-ibase+1;
      for( is=0; is<isize; is+=nblok1 ) {
        for( i=0; i<std::min(isize-is,nblok1); i++ ) {
          xig[iblok*nblok1+i] = xi[ibase+is+i];
          yig[iblok*nblok1+i] = yi[ibase+is+i];
          zig[iblok*nblok1+i] = zi[ibase+is+i];
          gxig[iblok*nblok1+i] = gxi[ibase+is+i];
          gyig[iblok*nblok1+i] = gyi[ibase+is+i];
          gzig[iblok*nblok1+i] = gzi[ibase+is+i];
          vig[iblok*nblok1+i] = 0;
        }
        for( i=isize-is; i<nblok1; i++ ) {
          xig[iblok*nblok1+i] = 0;
          yig[iblok*nblok1+i] = 0;
          zig[iblok*nblok1+i] = 0;
          gxig[iblok*nblok1+i] = 0;
          gyig[iblok*nblok1+i] = 0;
          gzig[iblok*nblok1+i] = 0;
          vig[iblok*nblok1+i] = 0;
        }
        nvecd[iblok*mblok+10] = nfi[ii];
        nvecd[iblok*mblok+11] = jbase;
        nvecd[iblok*mblok+12] = jsize;
        op += nblok1*jsize;
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
    xmind = xmin;
    ymind = ymin;
    zmind = zmin;
    l2pgpu_(nvecd,&op,&rb,&xmind,&ymind,&zmind,
            xig,yig,zig,gxig,gyig,gzig,vig,
            arex,aimx,fac);
    iblok = 0;
    for( ii=istagp[icall]; ii<=iendgp[icall]; ii++ ) {
      ibase = ndi[0][ii];
      isize = ndi[1][ii]-ibase+1;
      for( is=0; is<isize; is+=nblok1 ) {
        for( i=0; i<std::min(isize-is,nblok1); i++ ) {
          vi[ibase+is+i] += vig[iblok*nblok1+i];
        }
        iblok++;
      }
    }
  }

  toc = get_time();
  tfmm[22] += toc-tic;
}
