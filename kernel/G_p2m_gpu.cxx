#include "../include/constants.h"

extern float *xj,*yj,*zj,*gxj,*gyj,*gzj,*vj;
extern int **ndj,*nfj,*nc;
extern std::complex<double> (*bx)[mpsym],*bnm,*bth;
extern float *fac;
extern double *tfmm;
extern int *istagp,*iendgp,*nvecd;
extern float *xjg,*yjg,*zjg,*gxjg,*gyjg,*gzjg,*vjg;
extern float *brex,*bimx;

extern void p2mgpu_(int*, double*, double*, double*, double*, double*,
	float*, float*, float*, float*, float*, float*, float*,
	float*, float*, float*);

void G_p2m(int nmp, int mp, int lbj) {
  int kl,idev,mblok,ni,nj,ncall,jj,icall,iblok,jc,jbase,j,jsize,jm;
  double rb,op,xmind,ymind,zmind,tic2;
  tic = get_time();

  tic2 = get_time();
  kl = 1 << lmax;
  rb = rd/kl;
  idev = myrank%ngpu;
  mblok = 3;

  ni = 0;
  nj = 0;
  ncall = 0;
  istagp[0] = 0;
  for( jj=0; jj<lbj; jj++ ) {
    ni += ((nmp+nblok1)/nblok1+1)*nblok1;
    nj += ndj[1][jj]-ndj[0][jj]+1;
    if( ni > nimax || nj > njmax ) {
      iendgp[ncall] = jj-1;
      ncall++;
      istagp[ncall] = jj;
      ni = ((nmp+nblok1)/nblok1+1)*nblok1;
      nj = ndj[1][jj]-ndj[0][jj]+1;
    }
  }
  iendgp[ncall] = lbj-1;
  if(lbj != 0) ncall++;
  toc = get_time();
  tfmm[28] += toc-tic2;

  for( icall=0; icall<ncall; icall++ ) {
    tic2 = get_time();
    iblok = 0;
    jc = 0;
    op = 0;
    for( jj=istagp[icall]; jj<=iendgp[icall]; jj++ ) {
      jbase = jc;
      for( j=ndj[0][jj]; j<=ndj[1][jj]; j++ ) {
        xjg[jc] = xj[j];
        yjg[jc] = yj[j];
        zjg[jc] = zj[j];
        gxjg[jc] = gxj[j];
        gyjg[jc] = gyj[j];
        gzjg[jc] = gzj[j];
        vjg[jc] = vj[j];
        jc++;
      }
      jsize = jc-jbase;
      nvecd[iblok*mblok+10] = nfj[jj];
      nvecd[iblok*mblok+11] = jbase;
      nvecd[iblok*mblok+12] = jsize;
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
    xmind = xmin;
    ymind = ymin;
    zmind = zmin;
    toc = get_time();
    tfmm[29] += toc-tic2;
    tic2 = get_time();
    p2mgpu_(nvecd,&op,&rb,&xmind,&ymind,&zmind,
            xjg,yjg,zjg,gxjg,gyjg,gzjg,vjg,
            brex,bimx,fac);
    toc = get_time();
    tfmm[30] += toc-tic2;
    tic2 = get_time();
    iblok = 0;
    for( jj=istagp[icall]; jj<=iendgp[icall]; jj++ ) {
      for( j=0; j<nmp; j++ ) {
        jm = iblok*nblok1+j;
        bx[jj][j] = std::complex<double>(brex[jm],bimx[jm]);
      }
      iblok++;
    }
    toc = get_time();
    tfmm[29] += toc-tic2;
  }

  toc = get_time();
  tfmm[18] += toc-tic;
}
