#define NSCAL 6
#define MPMAX 100
#define NBLOK0 128
#define NBLOK1 64
#define ROUNDUP0(n) ((((n)+NBLOK0-1)/NBLOK0+1)*NBLOK0)
#define ROUNDUP1(n) ((((n)+NBLOK1-1)/NBLOK1+1)*NBLOK1)

#define LARGE_SHIFT 21
#define LARGE (float) (3 << (LARGE_SHIFT-1))

int idev,jdev,iblok;
static int is_set=0;
static double tgpu[10];
unsigned int ms,mn,mi,mj,mk,ml;
static unsigned int ms_a=0,mn_a=0,mi_a=0,mj_a=0,mk_a=0,ml_a=0;
int *nvec;
float *ivec;
float *jvec;
float *kvec;
float *lvec;
float *scald;
static int *nveg[4];
static float *iveg[4];
static float *jveg[4];
static float *kveg[4];
static float *lveg[4];

__device__ __constant__ float scal[NSCAL];

typedef struct {
  float hs;
  float ls;
} SS;
