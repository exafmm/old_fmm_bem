#include "../include/constants.h"

extern float *sortd;
extern double *tfmm;

void sortvar(int n0, int n1, float* var, int* nvar) {
  int i;
  tic = get_time();
  for( i=n0; i<n1; i++ ) {
    sortd[i] = var[nvar[i]];
  }
  for( i=n0; i<n1; i++ ) {
    var[i] = sortd[i];
  }
  toc = get_time();
  tfmm[8] += toc-tic;
}

void unsortvar(int n0, int n1, float* var, int* nvar) {
  int i;
  tic = get_time();
  for( i=n0; i<n1; i++ ) {
    sortd[nvar[i]] = var[i];
  }
  for( i=n0; i<n1; i++ ) {
    var[i] = sortd[i];
  }
  toc = get_time();
  tfmm[25] += toc-tic;
}
