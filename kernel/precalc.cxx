#include "../include/constants.h"

extern float *sr,*fac,*anm;
extern std::complex<double> (*bx)[mpsym],*ynm,***dnm;
extern double *tfmm;

extern double factorial(int);
extern void boxc(int, int, int*);

void precalc(int& nmp, int mp) {
  int n,m,nm,nabsm,j,k,jk,nk,jkn,jnk,npn,nmn,npm,nmm,nmk,i,nmk1,nm1k,nmk2,nc[3];
  double anmk[2][mpmax*mpmax],dnmd[mpmax*mpmax],fnmm,fnpm,fnma,fnpa,fad,pn,p,p1,p2;
  double anmd,anmkd,xijc,yijc,zijc,rh,al,be,sc,ank,ek;
  std::complex<double> expbe[mpmax],eim(0,1),cnm;
  tic = get_time();

  nmp = mp*(mp+1)/2;

  for( n=0; n<2*mp; n++ ) {
    for( m=-n; m<=n; m++ ) {
      nm = n*n+n+m;
      nabsm = abs(m);
      fnmm = factorial(n-m);
      fnpm = factorial(n+m);
      fnma = factorial(n-nabsm);
      fnpa = factorial(n+nabsm);
      fac[nm] = sqrt(fnma/fnpa);
      fad = sqrt(fnmm*fnpm);
      anm[nm] = pow(-1.0,n)/fad;
    }
  }

  for( j=0; j<mp; j++) {
    for( k=-j; k<=j; k++ ){
      jk = j*j+j+k;
      for( n=abs(k); n<mp; n++ ) {
        nk = n*n+n+k;
        jkn = jk*mp*mp+nk;
        jnk = (j+n)*(j+n)+j+n;
        sr[jkn] = pow(-1.0,j+k)*anm[nk]*anm[jk]/anm[jnk];
      }
    }
  }

  pn = 1;
  for( m=0; m<2*mp; m++ ) {
    p = pn;
    npn = m*m+2*m;
    nmn = m*m;
    ynm[npn] = fac[npn]*p;
    ynm[nmn] = conj(ynm[npn]);
    p1 = p;
    p = (2*m+1)*p;
    for( n=m+1; n<2*mp; n++ ) {
      npm = n*n+n+m;
      nmm = n*n+n-m;
      ynm[npm] = fac[npm]*p;
      ynm[nmm] = conj(ynm[npm]);
      p2 = p1;
      p1 = p;
      p = ((2*n+1)*p1-(n+m)*p2)/(n-m+1);
    }
    pn = 0;
  }

  for( n=0; n<mp; n++ ) {
    for( m=1; m<=n; m++ ) {
      anmd = n*(n+1)-m*(m-1);
      for( k=1-m; k<m; k++ ) {
        nmk = (4*n*n*n+6*n*n+5*n)/3+m*(2*n+1)+k;
        anmkd = ((double) (n*(n+1)-k*(k+1)))/(n*(n+1)-m*(m-1));
        anmk[0][nmk] = -(m+k)/sqrt(anmd);
        anmk[1][nmk] = sqrt(anmkd);
      }
    }
  }

  for( i=0; i<nrbm; i++ ) {
    boxc(i,3,nc);
    xijc = nc[0]-3;
    yijc = nc[1]-3;
    zijc = nc[2]-3;
    rh = sqrt(xijc*xijc+yijc*yijc+zijc*zijc)+eps;
    al = acos(zijc/rh);
    if( std::abs(xijc)+std::abs(yijc) < eps ) {
      be = 0;
    } else if( std::abs(xijc) < eps ) {
      be = yijc/std::abs(yijc)*pi*0.5;
    } else if( xijc > 0 ) {
      be = atan(yijc/xijc);
    } else {
      be = atan(yijc/xijc)+pi;
    }

    sc = sin(al)/(1+cos(al));
    for( n=0; n<4*mp-3; n++ ) {
      expbe[n] = exp((n-2*mp+2)*be*eim);
    }

    for( n=0; n<mp; n++ )  {
      nmk = (4*n*n*n+6*n*n+5*n)/3+n*(2*n+1)+n;
      dnmd[nmk] = pow(cos(al*0.5),2*n);
      for( k=n; k>=1-n; k-- ) {
        nmk = (4*n*n*n+6*n*n+5*n)/3+n*(2*n+1)+k;
        nmk1 = (4*n*n*n+6*n*n+5*n)/3+n*(2*n+1)+k-1;
        ank = ((double) n+k)/(n-k+1);
        dnmd[nmk1] = sqrt(ank)*tan(al*0.5)*dnmd[nmk];
      }
      for( m=n; m>=1; m-- ) {
        for( k=m-1; k>=1-m; k-- ){
          nmk = (4*n*n*n+6*n*n+5*n)/3+m*(2*n+1)+k;
          nmk1 = (4*n*n*n+6*n*n+5*n)/3+m*(2*n+1)+k+1;
          nm1k = (4*n*n*n+6*n*n+5*n)/3+(m-1)*(2*n+1)+k;
          dnmd[nm1k] = anmk[1][nmk]*dnmd[nmk1]+anmk[0][nmk]*sc*dnmd[nmk];
        }
      }
    }

    for( n=1; n<mp; n++ ) {
      for( m=0; m<=n; m++ ) {
        for( k=-m; k<=-1; k++ ) {
          ek = pow(-1.0,k);
          nmk = (4*n*n*n+6*n*n+5*n)/3+m*(2*n+1)+k;
          nmk1 = (4*n*n*n+6*n*n+5*n)/3-k*(2*n+1)-m;
          dnmd[nmk] = ek*dnmd[nmk];
          dnmd[nmk1] = pow(-1.0,m+k)*dnmd[nmk];
        }
        for( k=0; k<=m; k++ ) {
          nmk = (4*n*n*n+6*n*n+5*n)/3+m*(2*n+1)+k;
          nmk1 = (4*n*n*n+6*n*n+5*n)/3+k*(2*n+1)+m;
          nmk2 = (4*n*n*n+6*n*n+5*n)/3-k*(2*n+1)-m;
          dnmd[nmk1] = pow(-1.0,m+k)*dnmd[nmk];
          dnmd[nmk2] = dnmd[nmk1];
        }
      }
    }

    for( n=0; n<mp; n++ ) {
      for( m=0; m<=n; m++ ) {
        for( k=-n; k<=n; k++ ) {
          nmk = (4*n*n*n+6*n*n+5*n)/3+m*(2*n+1)+k;
          nk = n*(n+1)+k;
          dnm[i][m][nk] = dnmd[nmk]*expbe[k+m+2*mp-2];
        }
      }
    }

    al = -al;
    be = -be;

    sc = sin(al)/(1+cos(al));
    for( n=0; n<4*mp-3; n++ ) {
      expbe[n] = exp((n-2*mp+2)*be*eim);
    }

    for( n=0; n<mp; n++ ) {
      nmk = (4*n*n*n+6*n*n+5*n)/3+n*(2*n+1)+n;
      dnmd[nmk] = pow(cos(al*0.5),2*n);
      for( k=n; k>=1-n; k-- ) {
        nmk = (4*n*n*n+6*n*n+5*n)/3+n*(2*n+1)+k;
        nmk1 = (4*n*n*n+6*n*n+5*n)/3+n*(2*n+1)+k-1;
        ank = ((double) n+k)/(n-k+1);
        dnmd[nmk1] = sqrt(ank)*tan(al*0.5)*dnmd[nmk];
      }
      for( m=n; m>=1; m-- ) {
        for( k=m-1; k>=1-m; k-- ) {
          nmk = (4*n*n*n+6*n*n+5*n)/3+m*(2*n+1)+k;
          nmk1 = (4*n*n*n+6*n*n+5*n)/3+m*(2*n+1)+k+1;
          nm1k = (4*n*n*n+6*n*n+5*n)/3+(m-1)*(2*n+1)+k;
          dnmd[nm1k] = anmk[1][nmk]*dnmd[nmk1]+anmk[0][nmk]*sc*dnmd[nmk];
        }
      }
    }

    for( n=1; n<mp; n++ ) {
      for( m=0; m<=n; m++ ) {
        for( k=-m; k<=-1; k++ ) {
          ek = pow(-1.0,k);
          nmk = (4*n*n*n+6*n*n+5*n)/3+m*(2*n+1)+k;
          nmk1 = (4*n*n*n+6*n*n+5*n)/3-k*(2*n+1)-m;
          dnmd[nmk] = ek*dnmd[nmk];
          dnmd[nmk1] = pow(-1.0,m+k)*dnmd[nmk];
        }
        for( k=0; k<=m; k++ ) {
          nmk = (4*n*n*n+6*n*n+5*n)/3+m*(2*n+1)+k;
          nmk1 = (4*n*n*n+6*n*n+5*n)/3+k*(2*n+1)+m;
          nmk2 = (4*n*n*n+6*n*n+5*n)/3-k*(2*n+1)-m;
          dnmd[nmk1] = pow(-1.0,m+k)*dnmd[nmk];
          dnmd[nmk2] = dnmd[nmk1];
        }
      }
    }

    for( n=0; n<mp; n++ ) {
      for( m=0; m<=n; m++ ) {
        for( k=-n; k<=n; k++ ) {
          nmk = (4*n*n*n+6*n*n+5*n)/3+m*(2*n+1)+k;
          nk = n*(n+1)+k;
          dnm[i+nrbm][m][nk] = dnmd[nmk]*expbe[k+m+2*mp-2];
        }
      }
    }
  }

  for( j=0; j<nbnet; j++ ) {
    for( i=0; i<nmp; i++ ) {
      bx[j][i] = 0;
    }
  }

  toc = get_time();
  tfmm[3] += toc-tic;
}
