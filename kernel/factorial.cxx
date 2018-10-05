double factorial(int n) {
  int i;
  double fac = 1;
  if( n > 0 ) {
    for( i=1; i<=n; i++ ) {
      fac = fac*i;
    }
  }
  return fac;
}
