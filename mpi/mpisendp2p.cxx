#include "mpi.h"
#include "../include/constants.h"

extern float *xj,*yj,*zj,*gxj,*gyj,*gzj,*vj;
extern int **ndj,**nsij,*kscnt,*lsdsp,*lscnt,*lrdsp,*lrcnt;
extern float *fsend,*frecv;
extern double *tfmm;

extern void mpialltoallvf(float*, int*, int*, float*, int*, int*, int);

void mpisendp2p(int mj, int np, float* var) {
  int ic,ii,i,j,k;
  tic = get_time();

  ic = 0;
  for( ii=0; ii<nprocs; ii++ ) {
    for( i=0; i<kscnt[ii]; i++ ) {
      j = nsij[ii][i];
      for( k=ndj[0][j]; k<=ndj[1][j]; k++ ) {
        fsend[ic] = var[k];
        ic++;
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  mpialltoallvf(fsend,lscnt,lsdsp,frecv,lrcnt,lrdsp,npmax);
  ic = mj;
  for( i=0; i<np; i++ ) {
    var[ic] = frecv[i];
    ic++;
  }

  toc = get_time();
  tfmm[15] += toc-tic;
}
