#include <iostream>

extern int ncheck;
extern int *na,*nb;

void sort(int np) {
  int msta[10000],mend[10000],nstack,ins,ista,iend,j,ma,mb,i,imed,il,ir,naa,nbb;

  nstack = 1;
  ins = 15;
  ista = 1;
  iend = np;
  while( nstack != 0 ) {
    if( iend-ista < ins ) {
      for( j=ista+1; j<=iend; j++ ) {
        ma = na[j-1];
        mb = nb[j-1];
        for( i=j-1; i>=ista; i-- ) {
          if( na[i-1] <= ma ) break;
          na[i] = na[i-1];
          nb[i] = nb[i-1];
        }
        na[i] = ma;
        nb[i] = mb;
      }
      ista = msta[nstack-1];
      iend = mend[nstack-1];
      nstack--;
    } else {
      imed = (ista+iend)/2;
      ma = na[imed-1];
      il = ista-1;
      ir = iend+1;
      while( 1 ) {
        while ( 1 ) {
          il++;
          if( na[il-1] >= ma ) break;
        }
        while ( 1 ) {
          ir--;
          if( na[ir-1] <= ma ) break;
        }
        if( ir <= il ) break;
        naa = na[il-1];
        nbb = nb[il-1];
        na[il-1] = na[ir-1];
        nb[il-1] = nb[ir-1];
        na[ir-1] = naa;
        nb[ir-1] = nbb;
      }
      nstack++;
      if( iend-il+1 >= ir-ista+1 ) {
        msta[nstack-1] = il;
        mend[nstack-1] = iend;
        iend = ir;
      } else {
        msta[nstack-1] = ista;
        mend[nstack-1] = ir;
        ista = il;
      }
    }
  }
}
