#include "../include/constants.h"

void mpirange(int n0, int n1, int& ista, int& iend) {
  int iwork1,iwork2;

  iwork1 = (n1-n0+1)/nprocs;
  iwork2 = (n1-n0+1)%nprocs;
  ista = myrank*iwork1+n0+std::min(myrank,iwork2);
  iend = ista+iwork1-1;
  if( iwork2 > myrank ) iend++;
}
