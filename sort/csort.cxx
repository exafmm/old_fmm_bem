#include "../include/constants.h"

extern int *na,*nb,*nc,*nd;

void sort(int np) {
  int i;

  for( i=0; i<nbmax; i++ ) nc[i] = 0;
  for( i=0; i<np; i++ ) nc[na[i]]++;
  for( i=1; i<nbmax; i++ ) nc[i] += nc[i-1];
  for( i=np-1; i>=0; i-- ) {
    nc[na[i]]--;
    nd[nc[na[i]]] = na[i];
    nb[nc[na[i]]] = i;
  }
  for( i=0; i<np; i++ ) na[i] = nd[i];
}
