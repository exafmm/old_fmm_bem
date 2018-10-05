#include "../include/constants.h"

extern double *tfmm;

extern void ijbox(int, int, int, int);

void ifbox(int mi, int mj, int nmp, int lbi, int lbj, int& lbjr, int lev, int neib) {
  ijbox(lbi,lbj,lev,neib);
  tic = get_time();
  lbjr = lbj;
  toc = get_time();
  tfmm[12] += toc-tic;
}
