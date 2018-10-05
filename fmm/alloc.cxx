#include "../include/constants.h"

extern float *gxj,*gyj,*gzj;
extern int *nbi,**ndi,*nfi,*nlbi;
extern int *nbj,**ndj,*nej,*nfj,*nlbj,*nfo;
extern int *nij,*njb,(*neij)[nebm],**nsij,*nfrom;
extern double *tfmm;
extern int *nek,*irank;
extern std::complex<double> (*ax)[mpsym],(*axo)[mpsym],(*bx)[mpsym];
extern std::complex<double> *bnm,*bth,*ynm,***dnm;
extern float *fac,*sr,*anm;
extern int *iexp1,*iexp2;
extern float (*ytop)[mp+1];
extern float (*m2e)[mp+1][mp+1];
extern double (*rdpi2)[mp+1][2*mp+1];
extern double (*rdmpi2)[mp+1][2*mp+1];
extern double (*rdsq3)[mp+1][2*mp+1];
extern double (*rdmsq3)[mp+1][2*mp+1];
extern int *ibase,*isize,*jbase,*jsize;
extern int *istagp,*iendgp,**jstagp,**jendgp,*njcall,*nvecd,*njj;
extern float *xig,*yig,*zig,*gxig,*gyig,*gzig,*vig;
extern float *xjg,*yjg,*zjg,*gxjg,*gyjg,*gzjg,*vjg;
extern float *arex,*aimx,*brex,*bimx;
extern int *nbuf,*nbufd;
extern float *fbuf,*fbufd;
extern double *dbuf,*dbufd;
extern std::complex<float> *cbuf,*cbufd;

extern void memoryuse();
extern void memoryfree();

void alloc(int ialloc) {
  int i,j,nebmp;
  tic = get_time();

  if( ialloc == 1 ) {
    nebmp = std::max(nebm,nprocs);
  
    ndi = new int* [2];
    for( i=0; i<2; i++ ) ndi[i] = new int [nbne];
    nfi = new int [nbnet];
    nlbi = new int [lmax];
    ndj = new int* [2];
    for( i=0; i<2; i++ ) ndj[i] = new int [nbne];
    nej = new int [nbmax];
    nfj = new int [nbnet];
    nlbj = new int [lmax];
    nfo = new int [nbnet];
    nij = new int [nbne];
    njb = new int [nbne];
    neij = new int [nbne][nebm];
    nsij = new int* [nebmp];
    for( i=0; i<nebmp; i++ ) nsij[i] = new int [nbne];
    nfrom = new int [nbne];
    mem = nbmax*4+nbne*7*4+nbnet*3*4+nebm*nbne*4+nebmp*nbne*4;
    memoryuse();
  
    ax = new std::complex<double> [nbne][mpsym];
    axo = new std::complex<double> [nbne][mpsym];
    bx = new std::complex<double> [nbnet][mpsym];
    bnm = new std::complex<double> [4*mpmax];
    bth = new std::complex<double> [mpmax];
    sr = new float [mpmax*mpmax];
    iexp1 = new int [nbnet];
    iexp2 = new int [nbnet];
    ytop = new float [mp+1][mp+1];
    m2e = new float [mp][mp+1][mp+1];
    rdpi2 = new double [mp+1][mp+1][2*mp+1];
    rdmpi2 = new double [mp+1][mp+1][2*mp+1];
    rdsq3 = new double [mp+1][mp+1][2*mp+1];
    rdmsq3 = new double [mp+1][mp+1][2*mp+1];
    ynm = new std::complex<double> [4*mpmax];
    dnm = new std::complex<double>** [2*nrbm];
    for( i=0; i<2*nrbm; i++ ) {
      dnm[i] = new std::complex<double>* [mpcmp];
      for( j=0; j<mpcmp; j++ ) dnm[i][j] = new std::complex<double> [mpmax];
    }
    fac = new float [4*mpmax];
    anm = new float [4*mpmax];

    mem = mpsym*nbne*2*16+mpsym*nbnet*16+mpmax*9*16+mpmax*mpmax*4+mp*(mp+1)*(mp+1)*4+mpmax*mpcmp*2*nrbm*16;
    memoryuse();
    m2lrp2p = 200;
  }
  if( ngpu != 0 ) {
    jbase = new int [nbne];
    jsize = new int [nbne];
    istagp = new int [nbne];
    iendgp = new int [nbne];
    jstagp = new int* [nebm];
    for( i=0; i<nebm; i++ ) jstagp[i] = new int [nbne];
    jendgp = new int* [nebm];
    for( i=0; i<nebm; i++ ) jendgp[i] = new int [nbne];
    njcall = new int [nbne];
    nvecd = new int [npmax];
    njj = new int [nbne];
    xig = new float [nimax];
    yig = new float [nimax];
    zig = new float [nimax];
    gxig = new float [nimax];
    gyig = new float [nimax];
    gzig = new float [nimax];
    vig = new float [nimax];
    xjg = new float [njmax];
    yjg = new float [njmax];
    zjg = new float [njmax];
    gxjg = new float [njmax];
    gyjg = new float [njmax];
    gzjg = new float [njmax];
    vjg = new float [njmax];
    arex = new float [nimax];
    aimx = new float [nimax];
    brex = new float [njmax];
    bimx = new float [njmax];
    mem = npmax*4+nbne*6*4+nebm*nbne*2*4+nimax*9*4+njmax*9*4;
    memoryuse();
    m2lrp2p = 5000;
  }

  toc = get_time();
  tfmm[2] += toc-tic;
}

void dealloc(int ialloc) {
  int i,j,nebmp;
  tic = get_time();

  if( ialloc == 1 ) {
    nebmp = std::max(nebm,nprocs);
  
    for( i=0; i<2; i++ ) delete[] ndi[i];
    delete[] ndi;
    delete[] nfi;
    delete[] nlbi;
    for( i=0; i<2; i++ ) delete[] ndj[i];
    delete[] ndj;
    delete[] nej;
    delete[] nfj;
    delete[] nlbj;
    delete[] nfo;
    delete[] nij;
    delete[] njb;
    delete[] neij;
    for( i=0; i<nebmp; i++ ) delete[] nsij[i];
    delete[] nsij;
    delete[] nfrom;
    delete[] nek;
    delete[] irank;
    mem = nbmax*3*4+nbne*7*4+nbnet*3*4+nebm*nbne*4+nebmp*nbne*4;
    memoryfree();

    delete[] ax;
    delete[] axo;
    delete[] bx;
    delete[] bnm;
    delete[] bth;
    delete[] sr;
    delete[] iexp1;
    delete[] iexp2;
    delete[] ytop;
    delete[] m2e;
    delete[] rdpi2;
    delete[] rdmpi2;
    delete[] rdsq3;
    delete[] rdmsq3;
    delete[] ynm;
    for( i=0; i<2*nrbm; i++ ) {
      for( j=0; j<mpcmp; j++ ) delete[] dnm[i][j];
      delete[] dnm[i];
    }
    delete[] dnm;
    delete[] fac;
    delete[] anm;
    mem = mpsym*nbne*2*16+mpsym*nbnet*16+mpmax*9*16+mpmax*mpmax*4+mp*(mp+1)*(mp+1)*4+mpmax*mpcmp*2*nrbm*16;
    memoryfree();
  }

  if( ngpu != 0 ) {
    delete[] jbase;
    delete[] jsize;
    delete[] istagp;
    delete[] iendgp;
    for( i=0; i<nebm; i++ ) delete[] jstagp[i];
    delete[] jstagp;
    for( i=0; i<nebm; i++ ) delete[] jendgp[i];
    delete[] jendgp;
    delete[] njcall;
    delete[] nvecd;
    delete[] njj;
    delete[] xig;
    delete[] yig;
    delete[] zig;
    delete[] gxig;
    delete[] gyig;
    delete[] gzig;
    delete[] vig;
    delete[] xjg;
    delete[] yjg;
    delete[] zjg;
    delete[] gxjg;
    delete[] gyjg;
    delete[] gzjg;
    delete[] vjg;
    delete[] arex;
    delete[] aimx;
    delete[] brex;
    delete[] bimx;
    mem = npmax*4+nbne*6*4+nebm*nbne*2*4+nimax*9*4+njmax*9*4;
    memoryfree();
  }

  toc = get_time();
  tfmm[27] += toc-tic;
}
