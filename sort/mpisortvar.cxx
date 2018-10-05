#include "mpi.h"
#include "../include/constants.h"

extern float *sortd;
extern double *tfmm;

extern void mpialltoallvf(float*, int*, int*, float*, int*, int*, int);

void mpisortvar(int* nsdsp, int* nscnt, int* nrdsp, int*nrcnt, float* var) {
  int n,i;
  tic = get_time();

  n = 0;
  for( i=0; i<nprocs; i++ ) n += nrcnt[i];
  mpialltoallvf(var,nscnt,nsdsp,sortd,nrcnt,nrdsp,npmax);
  for( i=0; i<n; i++ ) var[i] = sortd[i];

  toc = get_time();
  tfmm[9] += toc-tic;
}

void mpiunsortvar(int* nsdsp, int* nscnt, int* nrdsp, int*nrcnt, float* var) {
  int n,i;
  tic = get_time();

  n = 0;
  for( i=0; i<nprocs; i++ ) n += nscnt[i];
  mpialltoallvf(var,nrcnt,nrdsp,sortd,nscnt,nsdsp,npmax);
  for( i=0; i<n; i++ ) var[i] = sortd[i];

  toc = get_time();
  tfmm[26] += toc-tic;
}
