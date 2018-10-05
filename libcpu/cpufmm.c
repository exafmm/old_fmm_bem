#include <assert.h>
#include <math.h>
#include <stdlib.h>

// summation
#define AX      "%xmm0"
#define AY      "%xmm1"
#define AZ      "%xmm2"
#define PHI     "%xmm3"
// j particle
#define XJ      "%xmm4"
#define YJ      "%xmm5"
#define ZJ      "%xmm6"
#define MJ      "%xmm7"
// temporary
#define RINV    "%xmm8"
#define X2      "%xmm9"
#define Y2      "%xmm10"
#define Z2      "%xmm11"
// fixed i particle
#define XI      "%xmm12"
#define YI      "%xmm13"
#define ZI      "%xmm14"
#define R2      "%xmm15"

#define XORPS(a, b) asm volatile("xorps "  a  ","  b );
#define LOADPS(mem, reg) asm volatile("movaps %0, %"reg::"m"(mem));
#define STORPS(reg, mem) asm volatile("movaps %"reg " , %0"::"m"(mem));
#define MOVAPS(src, dst) asm volatile("movaps " src "," dst);
#define MOVQ(src, dst) asm volatile("movq " src "," dst);
#define BCAST0(reg) asm volatile("shufps $0x00, " reg ","  reg);
#define BCAST1(reg) asm volatile("shufps $0x55, " reg ","  reg);
#define BCAST2(reg) asm volatile("shufps $0xaa, " reg ","  reg);
#define BCAST3(reg) asm volatile("shufps $0xff, " reg ","  reg);
#define MULPS(src, dst) asm volatile("mulps " src "," dst);
#define ADDPS(src, dst) asm volatile("addps " src ","  dst);
#define SUBPS(src, dst) asm volatile("subps "  src "," dst);
#define RSQRTPS(src, dst) asm volatile("rsqrtps " src "," dst);
#define MOVHLPS(src, dst) asm volatile("movhlps " src "," dst);
#define DEBUGPS(reg)

#define ALIGN16 __attribute__ ((aligned(16)))

#define NJMAX (1<<24)

typedef double v2df __attribute__ ((vector_size(16)));
typedef float  v4sf __attribute__ ((vector_size(16)));
typedef int    v4si __attribute__ ((vector_size(16)));
typedef short  v8hi __attribute__ ((vector_size(16)));

typedef struct iptdata{
  float x[4];
  float y[4];
  float z[4];
  float eps2[4]; // not used in this implementation
} Ipdata ALIGN16;
typedef struct fodata{
  float ax[4];
  float ay[4];
  float az[4];
  float phi[4];
} Fodata ALIGN16;
typedef struct jpdata{
  float x, y, z, m;
} Jpdata ALIGN16;

void GravityKernel(Ipdata *ipdata, Fodata *fodata, Jpdata *jpdata, int nj){
  int j;
  assert(((unsigned long)jpdata & 15) == 0);
  assert(((unsigned long)ipdata & 15) == 0);
  assert(((unsigned long)fodata & 15) == 0);

  XORPS(AX, AX);               // AX = 0
  XORPS(AY, AY);               // AY = 0
  XORPS(AZ, AZ);               // AZ = 0
  XORPS(PHI, PHI);             // PHI = 0

  LOADPS(*ipdata->x, XI);      // XI = *ipdata->x
  LOADPS(*ipdata->y, YI);      // YI = *ipdata->y
  LOADPS(*ipdata->z, ZI);      // ZI = *ipdata->z
  LOADPS(*ipdata->eps2, R2);   // R2 = *ipdata->eps2

  LOADPS(jpdata[0], MJ);       // MJ = *jpdata->x,y,z,m
  MOVAPS(MJ, X2);              // X2 = MJ
  MOVAPS(MJ, Y2);              // Y2 = MJ
  MOVAPS(MJ, Z2);              // Z2 = MJ

  BCAST0(X2);                  // X2 = *jpdata->x
  BCAST1(Y2);                  // Y2 = *jpdata->y
  BCAST2(Z2);                  // Z2 = *jpdata->z
  BCAST3(MJ);                  // MJ = *jpdata->m

  SUBPS(XI, X2);               // X2 = X2 - XI
  SUBPS(YI, Y2);               // Y2 = Y2 - YI
  SUBPS(ZI, Z2);               // Z2 = Z2 - ZI

  MOVAPS(X2, XJ);              // XJ = X2
  MOVAPS(Y2, YJ);              // YJ = Y2
  MOVAPS(Z2, ZJ);              // ZJ = Z2

  MULPS(X2, X2);               // X2 = X2 * X2
  MULPS(Y2, Y2);               // Y2 = Y2 * Y2
  MULPS(Z2, Z2);               // Z2 = Z2 * Z2

  ADDPS(X2, R2);               // R2 = R2 + X2
  ADDPS(Y2, R2);               // R2 = R2 + Y2
  ADDPS(Z2, R2);               // R2 = R2 + Z2

  LOADPS(jpdata[1], X2);       // X2 = *jpdata->x,y,z,m
  MOVAPS(X2, Y2);              // Y2 = X2
  MOVAPS(X2, Z2);              // Z2 = X2
  for(j=0;j<nj;j++){
    RSQRTPS(R2, RINV);         // RINV = rsqrt(R2)
    jpdata++;
    LOADPS(*ipdata->eps2, R2); // R2 = *ipdata->eps2
    BCAST0(X2);                // X2 = *jpdata->x
    BCAST1(Y2);                // Y2 = *jpdata->y
    BCAST2(Z2);                // Z2 = *jpdata->z
    SUBPS(XI, X2);             // X2 = X2 - XI
    SUBPS(YI, Y2);             // Y2 = Y2 - YI
    SUBPS(ZI, Z2);             // Z2 = Z2 - ZI

    MULPS(RINV, MJ);           // MJ = MJ * RINV
    SUBPS(MJ, PHI);            // PHI = PHI - MJ
    MULPS(RINV, RINV);         // RINV = RINV * RINV
    MULPS(MJ, RINV);           // RINV = MJ * RINV
    LOADPS(jpdata[0], MJ);     // MJ = *jpdata->x,y,z,m
    BCAST3(MJ);                // MJ = *jpdata->m

    MULPS(RINV, XJ);           // XJ = XJ * RINV
    ADDPS(XJ, AX);             // AX = AX + XJ
    MOVAPS(X2, XJ);            // XJ = X2
    MULPS(X2, X2);             // X2 = X2 * X2
    ADDPS(X2, R2);             // R2 = R2 + X2
    LOADPS(jpdata[1], X2);     // X2 = *jpdata->x,y,z,m

    MULPS(RINV, YJ);           // YJ = YJ * RINV
    ADDPS(YJ, AY);             // AY = AY + YJ
    MOVAPS(Y2, YJ);            // YJ = Y2
    MULPS(Y2, Y2);             // Y2 = Y2 * Y2
    ADDPS(Y2, R2);             // R2 = R2 + Y2
    MOVAPS(X2, Y2);            // Y2 = X2

    MULPS(RINV, ZJ);           // ZJ = ZJ * RINV
    ADDPS(ZJ, AZ);             // AZ = AZ + ZJ
    MOVAPS(Z2, ZJ);            // ZJ = Z2
    MULPS(Z2, Z2);             // Z2 = Z2 * Z2
    ADDPS(Z2, R2);             // R2 = R2 + Z2
    MOVAPS(X2, Z2);            // Z2 = X2
  }
  STORPS(AX, *fodata->ax);     // AX = *fodata->ax
  STORPS(AY, *fodata->ay);     // AY = *fodata->ay
  STORPS(AZ, *fodata->az);     // AZ = *fodata->az
  STORPS(PHI, *fodata->phi);   // PHI = *fodata->phi
}

static inline void v4sf_transpose(
    v4sf *d0, v4sf *d1, v4sf *d2, v4sf *d3,
    v4sf  s0, v4sf  s1, v4sf  s2, v4sf  s3)
{
  *d0 = __builtin_ia32_unpcklps(
        __builtin_ia32_unpcklps(s0, s2),
        __builtin_ia32_unpcklps(s1, s3));
  *d1 = __builtin_ia32_unpckhps(
        __builtin_ia32_unpcklps(s0, s2),
        __builtin_ia32_unpcklps(s1, s3));
  *d2 = __builtin_ia32_unpcklps(
        __builtin_ia32_unpckhps(s0, s2),
        __builtin_ia32_unpckhps(s1, s3));
  *d3 = __builtin_ia32_unpckhps(
        __builtin_ia32_unpckhps(s0, s2),
        __builtin_ia32_unpckhps(s1, s3));
}

static inline void v4sf_store_sp(v4sf vec, float *d0, float *d1, float *d2, float *d3){
  *d0 = __builtin_ia32_vec_ext_v4sf(vec, 0);
  *d1 = __builtin_ia32_vec_ext_v4sf(vec, 1);
  *d2 = __builtin_ia32_vec_ext_v4sf(vec, 2);
  *d3 = __builtin_ia32_vec_ext_v4sf(vec, 3);
}

void p2pcpu_(int nvecd[],double *op,
  float xi[],float yi[],float zi[],float gxi[],float gyi[],float gzi[],float vi[],
  float xj[],float yj[],float zj[],float gxj[],float gyj[],float gzj[],float vj[])
{
  int ii,icell=-1,ij,nij,jbase,jsize,nj=0,j,i,off;
  int iblok=nvecd[1],mblok=nvecd[2];
  Ipdata iptcl;
  Fodata fout;
  Jpdata *jptcl;
  jptcl = (Jpdata *) malloc(sizeof(Jpdata)*NJMAX);
  for( ii=0; ii<iblok; ii++ ){
    if( icell != nvecd[ii*mblok+10] ) {
      icell = nvecd[ii*mblok+10];
      nij = nvecd[ii*mblok+11];
      nj = 0;
      for( ij=0; ij<nij; ij++ ){
        jbase=nvecd[ii*mblok+2*ij+12];\
        jsize=nvecd[ii*mblok+2*ij+13];\
        for( j=jbase; j<jbase+jsize; j++ ){
          *(v4sf *)(jptcl+nj) = (v4sf) {xj[j],yj[j],zj[j],vj[j]};
          nj++;
        }
      }
    }
    off = 4*ii;
    for(i=0;i<4;i++){
      iptcl.x[i] = xi[off+i];
      iptcl.y[i] = yi[off+i];
      iptcl.z[i] = zi[off+i];
      iptcl.eps2[i] = 1e-12;
    }
    GravityKernel(&iptcl, &fout, jptcl, nj);
    v4sf ax = *(v4sf *)(fout.ax);
    v4sf ay = *(v4sf *)(fout.ay);
    v4sf az = *(v4sf *)(fout.az);
    v4sf phi= -*(v4sf *)(fout.phi);
    v4sf f0, f1, f2, f3;
    v4sf_transpose(&f0, &f1, &f2, &f3, ax, ay, az, phi);
    v4sf_store_sp(f0, &gxi[off+0], &gyi[off+0], &gzi[off+0], &vi[off+0]);
    v4sf_store_sp(f1, &gxi[off+1], &gyi[off+1], &gzi[off+1], &vi[off+1]);
    v4sf_store_sp(f2, &gxi[off+2], &gyi[off+2], &gzi[off+2], &vi[off+2]);
    v4sf_store_sp(f3, &gxi[off+3], &gyi[off+3], &gzi[off+3], &vi[off+3]);
    for(i=0;i<4;i++){
      gxi[off+i] *= 0.25/M_PI;
      gyi[off+i] *= 0.25/M_PI;
      gzi[off+i] *= 0.25/M_PI;
      vi[off+i] *= 0.25/M_PI;
    }
  }
  free(jptcl);
}
