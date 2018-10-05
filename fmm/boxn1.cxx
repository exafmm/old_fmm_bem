#include <cmath>

void boxn1(int* nc, int& n, int lev) {
  int i,nfx,nfy,nfz;
  n = 0;
  for( i=0; i<lev; i++ ) {
    nfx = nc[0]%2;
    nc[0] /= 2;
    n += nfx*(1 << (3*i+1));

    nfy = nc[1]%2;
    nc[1] /= 2;
    n += nfy*(1 << (3*i));

    nfz = nc[2]%2;
    nc[2] /= 2;
    n += nfz*(1 << (3*i+2));
  }
}
