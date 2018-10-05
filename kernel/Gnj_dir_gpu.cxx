#include "../include/constants.h"

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi,*vi;
extern float *xj,*yj,*zj,*gxj,*gyj,*gzj,*vj;
extern int *nvecd;
extern float *xig,*yig,*zig,*gxig,*gyig,*gzig,*vig;
extern float *xjg,*yjg,*zjg,*gxjg,*gyjg,*gzjg,*vjg;

extern void p2pgpu_(int*, double*, float*, float*, float*, float*, float*, float*, float*,
        float*, float*, float*, float*, float*, float*, float*);

void Gnj_dir(int n0, int n1, int n2, int n3) {
  int idev,i,nicall,njcall,icall,iwork1,iwork2,ista,iend,ibase,isize,iblok,is,mblok;
  int jcall,jwork1,jwork2,jsta,jend,jbase,jsize,nj=0;
  double op;

  idev = myrank%ngpu;
  for( i=n0; i<=n1; i++ ) vi[i] = 0;
  nicall = (n1-n0+1)/nimax+1;
  njcall = (n3-n2+1)/njmax+1;
  if( n0>n1 || n2>n3 ) nicall = 0;
  for( icall=0; icall<nicall; icall++ ) {
    iwork1 = (n1-n0+1)/nicall;
    iwork2 = (n1-n0+1)%nicall;
    ista = icall*iwork1+n0+std::min(icall,iwork2);
    iend = ista+iwork1-1;
    if( iwork2 > icall ) iend++;
    ibase = ista;
    isize = iend-ibase+1;
    iblok = 0;
    for( is=0; is<isize; is+=nblok0 ) {
      for( i=0; i<std::min(isize-is,nblok0); i++ ) {
        xig[iblok*nblok0+i] = xi[ibase+is+i];
        yig[iblok*nblok0+i] = yi[ibase+is+i];
        zig[iblok*nblok0+i] = zi[ibase+is+i];
        gxig[iblok*nblok0+i] = gxi[ibase+is+i];
        gyig[iblok*nblok0+i] = gyi[ibase+is+i];
        gzig[iblok*nblok0+i] = gzi[ibase+is+i];
        vig[iblok*nblok0+i] = vi[ibase+is+i];
      }
      for( i=isize-is; i<nblok0; i++ ) {
        xig[iblok*nblok0+i] = 0;
        yig[iblok*nblok0+i] = 0;
        zig[iblok*nblok0+i] = 0;
        gxig[iblok*nblok0+i] = 0;
        gyig[iblok*nblok0+i] = 0;
        gzig[iblok*nblok0+i] = 0;
        vig[iblok*nblok0+i] = 0;
      }
      iblok++;
    }
    mblok = 3;
    for( jcall=0; jcall<njcall; jcall++ ) {
      jwork1 = (n3-n2+1)/njcall;
      jwork2 = (n3-n2+1)%njcall;
      jsta = jcall*jwork1+n2+std::min(jcall,jwork2);
      jend = jsta+jwork1;
      if( jwork2 > jcall ) jend++;
      jbase = jsta;
      jsize = jend-jbase;
      for( i=0; i<iblok; i++ ) {
        nvecd[i*mblok+10] = 1;
        nvecd[i*mblok+11] = 0;
        nvecd[i*mblok+12] = jsize;
      }
      for( i=jsta; i<jend; i++ ) {
        nj = i-jsta;
        xjg[nj] = xj[i];
        yjg[nj] = yj[i];
        zjg[nj] = zj[i];
        gxjg[nj] = gxj[i];
        gyjg[nj] = gyj[i];
        gzjg[nj] = gzj[i];
        vjg[nj] = vj[i];
      }
      nj++;
      op = (double) isize*jsize;
      nvecd[0] = idev;
      nvecd[1] = iblok;
      nvecd[2] = mblok;
      nvecd[3] = nj;
      nvecd[4] = 2;
      nvecd[5] = myrank;
      p2pgpu_(nvecd,&op,xig,yig,zig,gxig,gyig,gzig,vig,
              xjg,yjg,zjg,gxjg,gyjg,gzjg,vjg);
      iblok = 0;
      for( is=0; is<isize; is+=nblok0 ) {
        for( i=0; i<std::min(isize-is,nblok0); i++ ) {
          vi[ibase+is+i] += vig[iblok*nblok0+i];
        }
        iblok++;
      }
    }
  }
}
