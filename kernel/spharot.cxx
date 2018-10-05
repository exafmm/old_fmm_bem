#include "../include/constants.h"

void spharot(std::complex<double>* yx, std::complex<double>* yxr, int mp, std::complex<double>** dnm) {
  int n,m,nms,k,nk,nks;
  std::complex<double> yxd;

  for( n=0; n<mp; n++ ) {
    for( m=0; m<=n; m++ ) {
      nms = n*(n+1)/2+m;
      yxd = 0;
      for( k=-n; k<=-1; k++ ) {
        nk = n*(n+1)+k;
        nks = n*(n+1)/2-k;
        yxd += dnm[m][nk]*conj(yx[nks]);
      }
      for( k=0; k<=n; k++ ) {
        nk = n*(n+1)+k;
        nks = n*(n+1)/2+k;
        yxd += dnm[m][nk]*yx[nks];
      }
      yxr[nms] = yxd;
    }
  }
}
