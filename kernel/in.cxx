#include "../include/constants.h"

void ribesl(double, double, int, int, double*, int);
double dgamma(double);

void in(double scal, double x, const int nb, double* b, int ncalc) {
  int i,ize;
  double cons, ensig, enmten, halfpi, xscal, term1, term2, alpha;
  ensig = 1.e-4;
  enmten = 1.e-300;

  if( x < 0.e0 ) printf("error in input for in.cxx\n");

  if( x <= ensig ) {
    xscal = x/scal;
    term1 = 1.0e0;
    term2 = 0.5e0*x*x;
    b[0] = term1*(1+term2/3);
    for( i=1; i<=nb; i++ ) {
      term1 *= xscal/double(2*i+1);
      if( term1 <= enmten ) term1 = 0;
      b[i] = term1*(1+term2/double(2*i+3));
    }
    ncalc = nb+1;
  } else if( x > 100 ) {
    for( i=0; i<nb+1; i++ ) b[i] = 0;
    ncalc = nb+1;
  } else {
    halfpi = atan(1.0)*2.0;
    cons = sqrt(halfpi/x);
    alpha = 0.5;
    ize = 1;
    ribesl(x,alpha,nb+1,ize,b,ncalc);
    for( i=0; i<nb+1; i++ ) {
      b[i] *= cons;
      cons /= scal;
      if( std::abs(b[i]) < enmten ) cons = 0;
    }
  }
}

void ribesl(double x, double alpha, int nb, int ize, double* b, int ncalc) {
  int k,l,magx,n,nbmx,nend,nsig,nstart;
  double cons,em,empal,emp2al,en,enmten,ensig,enten,exparg,half,halfx,
    one,p,plast,pold,psave,psavel,rtnsig,sums,tempa,tempb,tempc,test,tover,two,xlarge,zero;
  one = 1e0;
  two = 2e0;
  zero = 0e0;
  half = .5e0;
  cons = 1.585e0;
  nsig = 16;
  xlarge = 1e4;
  exparg = 709e0;
  enten = 1e308;
  ensig = 1e16;
  rtnsig = 1e-4;
  enmten = 8.9e-308;

  if( ( nb > 0 ) && ( x >= zero ) && ( alpha >= zero ) && ( alpha < one ) &&
    ( ( ( ize == 1 ) && ( x <= exparg ) ) || ( ( ize == 2 ) && ( x <= xlarge ) ) ) ) {
    ncalc = nb;
    magx = int(x);
    if( x >= rtnsig ) {
      nbmx = nb-magx;
      n = magx+1;
      en = double(n+n) + (alpha+alpha);
      plast = one;
      p = en/x;
      test = ensig + ensig;
      if( 2*magx > 5*nsig ) {
        test = sqrt(test*p);
      } else {
        test /= pow(cons,magx);
      }
      if( nbmx >= 3 ) {
        tover = enten/ensig;
        nstart = magx+2;
        nend = nb-1;
        for( k=nstart; k<=nend; k++ ) {
          n = k;
          en += two;
          pold = plast;
          plast = p;
          p = en*plast/x+pold;
          if( p > tover ) {
            tover = enten;
            p /= tover;
            plast /= tover;
            psave = p;
            psavel = plast;
            nstart = n+1;
l60:        n++;
            en += two;
            pold = plast;
            plast = p;
            p = en*plast/x+pold;
            if( p <= one ) goto l60;
            tempb = en/x;
            test = pold*plast/ensig;
            test *= half-half/(tempb*tempb);
            p = plast*tover;
            n--;
            en -= two;
            nend = std::min(nb,n);
            for( l=nstart; l<=nend; l++ ) {
              ncalc = l;
              pold = psavel;
              psavel = psave;
              psave = en*psavel/x+pold;
              if( psave*psavel > test ) goto l90;
            }
            ncalc = nend+1;
l90:        ncalc--;
            goto l120;
          }
        }
        n = nend;
        en = double(n+n)+(alpha+alpha);
        test = std::max(test,sqrt(plast*ensig)*sqrt(p+p));
      }
l110: n++;
      en += two;
      pold = plast;
      plast = p;
      p = en*plast/x+pold;
      if( p < test ) goto l110;
l120: n++;
      en += two;
      tempb = zero;
      tempa = one/p;
      em = double(n) - one;
      empal = em+alpha;
      emp2al = (em-one)+(alpha+alpha);
      sums = tempa*empal*emp2al/em;
      nend = n-nb;
      if( nend < 0 ) {
        b[n-1] = tempa;
        nend = -nend;
        for( l=1; l<=nend; l++ ) b[n+l-1] = zero;
      } else {
        if( nend > 0 ) {
          for( l=1; l<=nend; l++ ) {
            n--;
            en -= two;
            tempc = tempb;
            tempb = tempa;
            tempa = (en*tempb)/x+tempc;
            em -= one;
            emp2al -= one;
            if( n == 1 ) goto l150;
            if( n == 2 ) emp2al = one;
            empal -= one;
            sums = (sums+tempa*empal)*emp2al/em;
          }
        }
l150:   b[n-1] = tempa;
        if( nb <= 1 ) {
          sums = (sums+sums)+tempa;
          goto l230;
        }
        n--;
        en -= two;
        b[n-1] = (en*tempa)/x+tempb;
        if( n == 1 ) goto l220;
        em -= one;
        emp2al -= one;
        if( n == 2 ) emp2al = one;
        empal -= one;
        sums = (sums+b[n-1]*empal)*emp2al/em;
      }
      nend = n-2;
      if( nend > 0 ) {
        for( l=1; l<=nend; l++ ) {
          n--;
          en -= two;
          b[n-1] = (en*b[n])/x+b[n+1];
          em -= one;
          emp2al -= one;
          if( n == 2 ) emp2al = one;
          empal -= one;
          sums = (sums+b[n-1]*empal)*emp2al/em;
        }
      }
      b[0] = two*empal*b[1]/x+b[2];
l220: sums = (sums+sums)+b[0];
l230: if( alpha != zero ) sums *= dgamma(one+alpha)*pow(x*half,-alpha);
      if( ize == 1 ) sums *= exp(-x);
      tempa = enmten;
      if( sums > one ) tempa *= sums;
      for( n=1; n<=nb; n++ ) {
        if( b[n-1] < tempa ) b[n-1] = zero;
        b[n-1] /= sums;
      }
      return;
    } else {
      tempa = one;
      empal = one+alpha;
      halfx = zero;
      if( x > enmten ) halfx = half*x;
      if( alpha != zero ) tempa = pow(halfx,alpha)/dgamma(empal);
      if( ize == 2 ) tempa *= exp(-x);
      tempb = zero;
      if( ( x + one ) > one ) tempb = halfx*halfx;
      b[0] = tempa+tempa*tempb/empal;
      if( ( x != zero ) && ( b[0] == zero ) ) ncalc = 0;
      if( nb > 1 ) {
        if( x == zero ) {
          for( n=2; n<=nb; n++ ) b[n-1] = zero;
        } else {
          tempc = halfx;
          tover = (enmten+enmten)/x;
          if( tempb != zero ) tover = enmten/tempb;
          for( n=2; n<=nb; n++ ) {
            tempa /= empal;
            empal += one;
            tempa *= tempc;
            if( tempa <= tover*empal ) tempa = zero;
            b[n-1] = tempa+tempa*tempb/empal;
            if( ( b[n-1] == zero ) && ( ncalc > n ) ) ncalc = n-1;
          }
        }
      }
    }
  } else {
    ncalc = std::min(nb,0)-1;
  }
}

double dgamma(double x) {
  int i,n;
  bool parity;
  double eps,fact,half,one,pi,res,sqrtpi,sums,twelve,
    two,xbig,xden,xinf,xminin,xnum,y,y1,ysq,z,zero;
  double c[7]={-1.910444077728e-3,8.4171387781295e-4,-5.952379913043012e-4,7.93650793500350248e-4,
               -2.777777777777681622553e-3,8.333333333333333331554247e-2,5.7083835261e-3};
  double p[8]={-1.71618513886549492533811e0,2.47656508055759199108314e1,
               -3.79804256470945635097577e2,6.29331155312818442661052e2,
                8.66966202790413211295064e2,-3.14512729688483675254357e4,
               -3.61444134186911729807069e4,6.64561438202405440627855e4};
  double q[8]={-3.08402300119738975254353e1,3.15350626979604161529144e2,
               -1.01515636749021914166146e3,-3.10777167157231109440444e3,
                2.25381184209801510330112e4,4.75584627752788110767815e3,
               -1.34659959864969306392456e5,-1.15132259675553483497211e5};
  one = 1e0;
  half = .5e0;
  twelve = 12e0;
  two = 2e0;
  zero = 0e0;
  sqrtpi = 0.9189385332046727417803297e0;
  pi = 3.1415926535897932384626434e0;
  xbig = 171.624e0;
  xminin = 2.23e-308;
  eps = 2.22e-16;
  xinf = 1.79e308;
  parity = false;
  fact = one;
  n = 0;
  y = x;
  if( y <= zero ) {
    y = -x;
    y1 = floor(y);
    res = y-y1;
    if( res != zero ) {
      if( y1 != floor(y1*half)*two ) parity = true;
      fact = -pi/sin(pi*res);
      y += one;
    } else {
      res = xinf;
      goto l900;
    }
  }
  if( y < eps ) {
    if( y >= xminin ) {
      res = one/y;
    } else {
      res = xinf;
      goto l900;
    }
  } else if( y < twelve ) {
    y1 = y;
    if( y < one ) {
      z = y;
      y += one;
    } else {
      n = int(y)-1;
      y = y-double(n);
      z = y-one;
    }
    xnum = zero;
    xden = one;
    for( i=0; i<8; i++ ) {
      xnum = (xnum+p[i])*z;
      xden = xden*z+q[i];
    }
    res = xnum/xden+one;
    if( y1 < y ) {
      res /= y1;
    } else if( y1 > y ) {
      for( i=1; i<=n; i++ ) {
        res *= y;
        y += one;
      }
    }
  } else {
    if( y <= xbig ) {
      ysq = y*y;
      sums = c[6];
      for( i=0; i<6; i++ ) sums = sums/ysq+c[i];
      sums = sums/y-y+sqrtpi;
      sums = sums+(y-half)*log(y);
      res = exp(sums);
    } else {
      res = xinf;
      goto l900;
    }
  }
  if(parity) res = -res;
  if(fact != one ) res = fact/res;
l900: return res;
}
