#include <complex>

// coordinates, strengths

float *xi,*yi,*zi;
float *gxi,*gyi,*gzi,*vi;
float *xj,*yj,*zj;
float *gxj,*gyj,*gzj,*vj;
float *xq,*yq,*zq,*vq;
float *xo,*yo,*zo;
float *gxo,*gyo,*gzo;
float *uo,*vo,*wo;

// BEM quadratures and node values

int *face,*idx;
float *vert,(*xvert)[6],(*yvert)[6],(*zvert)[6];

float (*ptr)[3];
int (*ntr)[6],(*nne)[20],*im;
float (*vna)[3],(*vta)[3],(*vsa)[3];
float (*p1)[13],(*p2)[13],(*p3)[13];
float (*p4)[13],(*p5)[13],(*p6)[13];
float (*hsg)[13],*are,*crv;
double *si,*un;
float *xiq,*etq,*wq;

// FMM indexing

int *nbi,**ndi,*nfi,*nlbi;
int *nbj,**ndj,*nej,*nfj,*nlbj;
int *neo,*nfo,*nfn;
int *nij,*njb,(*neij)[nebm],**nsij;
int *ncnt,*nfrom;

// FMM moments and translation components

std::complex<double> (*ax)[mpsym],(*axo)[mpsym],(*bx)[mpsym];
std::complex<double> *bnm,*bth,*ynm,***dnm;
float *fac,*sr,*anm,*ank;

// FMM Helmholtz variables

int *iexp1,*iexp2;
float (*ytop)[mp+1];
float (*m2e)[mp+1][mp+1];
double (*rdpi2)[mp+1][2*mp+1];
double (*rdmpi2)[mp+1][2*mp+1];
double (*rdsq3)[mp+1][2*mp+1];
double (*rdmsq3)[mp+1][2*mp+1];
std::complex<double> *fexpe,*fexpo,*fexpb;
std::complex<double> (*xs)[3],(*ys)[3];
std::complex<double> **lexp1,**lexp2;
double (*zs)[3];

// FMM auxiliary variables

float *xd,*sortd;
int *nsortd;
double *tfmm,(*tall)[100];

// ParMETIS variables

int *nek,*irank;

// GPU indexing

int *ibase,*isize,*jbase,*jsize;
int *istagp,*iendgp,*jstagp,*jendgp;
int *njcall;
int *nvecd,*njj;

// GPU coordinates, strengths

float *xig,*yig,*zig;
float *gxig,*gyig,*gzig,*vig;
float *xjg,*yjg,*zjg;
float *gxjg,*gyjg,*gzjg,*vjg;

// GPU multipole & local expansions

float *arex,*aimx,*brex,*bimx;
float *ynmre,*ynmim;
float *dnmre,*dnmim;

// MPI indexing & variables

int *isdsp,*iscnt,*irdsp,*ircnt,*isort;
int *jsdsp,*jscnt,*jrdsp,*jrcnt,*jsort;
int *ksdsp,*kscnt,*krdsp,*krcnt;
int *lsdsp,*lscnt,*lrdsp,*lrcnt;
int *nsdsp,*nscnt,*nrdsp,*nrcnt;
int *nsend,*nrecv,*nbuf,*nbufd;
float *fsend,*frecv,*fbuf,*fbufd;
double *dsend,*drecv,*dbuf,*dbufd;
std::complex<float> *csend,*crecv,*cbuf,*cbufd;

// Sort variables

int *na,*nb,*nc,*nd;
