#include "../include/constants.h"

extern double *tfmm;

void nlevel(int nmax){
//  float level_switch[6]={1e4,4e4,3e5,2e6,1.5e7,1.3e8}; // cpu
//  float level_switch[6]={1e5,7e5,7e6,5e7,3e8,1.5e9}; // gpu
  float level_switch[6]={1e5,5e5,3e6,2e7,1.5e8,1.3e9}; // gpu p=5
  tic = get_time();

  lmax = 1;
  if( nmax < level_switch[0] ) {
    lmax += 1;
  } else if( nmax < level_switch[1] ) {
    lmax += 2;
  } else if( nmax < level_switch[2] ) {
    lmax += 3;
  } else if( nmax < level_switch[3] ) {
    lmax += 4;
  } else if( nmax < level_switch[4] ) {
    lmax += 5;
  } else if( nmax < level_switch[5] ) {
    lmax += 6;
  } else {
    lmax += 7;
  }
  printf("lmax   : %d\n",lmax);

  toc = get_time();
  tfmm[0] += toc-tic;
}
