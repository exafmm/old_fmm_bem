#include "mpi.h"
#include "../include/constants.h"

extern double *tfmm;

void nlevel(int mmax){
  int nmax;
//  float level_switch[7]={3e4,2e5,1.5e6,1.3e7,1e8,7e8,5e9}; // cpu 
//  float level_switch[7]={1e5,7e5,7e6,5e7,3e8,1.5e9,1.3e10}; // gpu
  float level_switch[7]={1e5,7e5,5e6,3e7,2e8,1.5e9,1.3e10}; // gpu p=5
//  float level_switch[7]={3e3,1.3e4,1e5,2e5,3e5,6e5,1e8}; // cpu bem
//  float level_switch[7]={2e4,1.3e5,3e5,5e5,8e5,1e7,2e8}; // gpu bem
  tic = get_time();

  MPI_Allreduce(&mmax,&nmax,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

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
  } else if( nmax < level_switch[6] ) {
    lmax += 7;
  } else {
    lmax += 8;
  }
  while( nprocs > ( 1 << 3*lmax ) ) lmax++;
  if( myrank == 0 ) printf("lmax   : %d\n",lmax);

  toc = get_time();
  tfmm[0] += toc-tic;
}
