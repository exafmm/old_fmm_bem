#include <cmath>

void boxc(int nf, int ndi, int* nc) {
  int i,nb,k,j,ndu;

  for( i=0; i<3; i++ ) nc[i] = 0;
  nb = nf;
  k = 0;
  i = 0;
  while( nb != 0 ) {
    j = ndi-k-1;
    nc[j] += (nb%2)*(1 << i);
    nb /= 2;
    k = (k+1)%ndi;
    if( k == 0 ) i++;
  }

  if( ndi == 3 ) {
    ndu = nc[0];
    nc[0] = nc[1];
    nc[1] = nc[2];
    nc[2] = ndu;
  }
}
