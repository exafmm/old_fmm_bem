#include "../include/constants.h"

extern float (*m2e)[mp+1][mp+1];
extern float (*ytop)[mp+1];
extern double (*rdpi2)[mp+1][2*mp+1];
extern double (*rdmpi2)[mp+1][2*mp+1];
extern double (*rdsq3)[mp+1][2*mp+1];
extern double (*rdmsq3)[mp+1][2*mp+1];
extern std::complex<double> *fexpe,*fexpo,*fexpb;
extern std::complex<double> (*xs)[3],(*ys)[3];
extern double (*zs)[3];

void fstrtn(double sqc[][4*mp+1], double d[][mp+1][2*mp+1], double theta) {
  int l,m,n,ij,im,imp,labs;
  double ctheta, stheta, hsthta;
  double cthtap, cthtan, precis=1e-19;
  double ww=0.7071067811865476;
  double facts[200];

  ctheta = cos(theta);
  if( fabs(ctheta) <= precis ) ctheta = 0;
  stheta = sin(-theta);
  if( fabs(stheta) <= precis ) stheta = 0;
  hsthta = ww*stheta;
  cthtap = ww*(1+ctheta);
  cthtan =-ww*(1-ctheta);

  d[0][0][mp] = 1;

  for( ij=1; ij<mp+1; ij++ ) {
    for( im=-ij; im<=-1; im++ ) {
      d[ij][0][im+mp] = -sqc[ij-im][2]*d[ij-1][0][im+mp+1];
      if( im > 1-ij ) {
        d[ij][0][im+mp] += sqc[ij+im][2]*d[ij-1][0][im+mp-1];
      }
      d[ij][0][im+mp]=d[ij][0][im+mp]*hsthta;
      if( im > -ij ) {
        d[ij][0][im+mp] += d[ij-1][0][im+mp]*ctheta*sqc[ij+im][1]*sqc[ij-im][1];
      }
      d[ij][0][im+mp] /= ij;
    }

    d[ij][0][mp] = d[ij-1][0][mp]*ctheta;
    if( ij > 1 ) {
      d[ij][0][mp] += hsthta*sqc[ij][2]*(d[ij-1][0][mp-1]+d[ij-1][0][mp+1])/ij;
    }

    for( im=1; im<=ij; im++ ) {
      d[ij][0][im+mp] = -sqc[ij+im][2]*d[ij-1][0][im+mp-1];
      if( im < ij-1 ) {
        d[ij][0][im+mp] += sqc[ij-im][2]*d[ij-1][0][im+mp+1];
      }
      d[ij][0][im+mp] *= hsthta;
      if( im < ij ) {
        d[ij][0][im+mp] += d[ij-1][0][im+mp]*ctheta*sqc[ij+im][1]*sqc[ij-im][1];
      }
      d[ij][0][im+mp] /= ij;
    }

    for( imp=1; imp<=ij; imp++ ) {
      for( im=-ij; im<=-1; im++ ) {
        d[ij][imp][im+mp] = d[ij-1][imp-1][im+mp+1]*cthtan*sqc[ij-im][2];
        if( im > 1-ij ) {
          d[ij][imp][im+mp] -= d[ij-1][imp-1][im+mp-1]*cthtap*sqc[ij+im][2];
        }
        if( im > -ij ) {
          d[ij][imp][im+mp] += d[ij-1][imp-1][im+mp]*stheta*sqc[ij+im][1]*sqc[ij-im][1];
        }
        d[ij][imp][im+mp] *= ww/sqc[ij+imp][2];
      }

      d[ij][imp][mp] = ij*stheta*d[ij-1][imp-1][mp];
      if( ij > 1 ) {
        d[ij][imp][mp] -= sqc[ij][2]*(d[ij-1][imp-1][mp-1]*cthtap+d[ij-1][imp-1][mp+1]*cthtan);
      }
      d[ij][imp][mp] *= ww/sqc[ij+imp][2];

      for( im=1; im<=ij; im++ ) {
        d[ij][imp][im+mp] = d[ij-1][imp-1][im+mp-1]*cthtap*sqc[ij+im][2];
        if( im < ij-1 ) {
          d[ij][imp][im+mp] -= d[ij-1][imp-1][im+mp+1]*cthtan*sqc[ij-im][2];
        }
        if( im < ij ) {
          d[ij][imp][im+mp] += d[ij-1][imp-1][im+mp]*stheta*sqc[ij+im][1]*sqc[ij-im][1];
        }
        d[ij][imp][im+mp] *= ww/sqc[ij+imp][2];
      }
    }
  }
  facts[0] = 1;
  for( n=1; n<2*mp+1; n++ ) {
    facts[n] = facts[n-1]*n;
  }
  for( n=0; n<mp+1; n++ ) {
    for( m=0; m<=n; m++ ) {
      for( l=-n; l<=n; l++ ) {
        labs = abs(l);
        d[n][m][mp+l] *= sqrt(facts[n+m]/facts[n+labs]*facts[n-labs]/facts[n-m]);
      }
    }
  }
}

void precalch(int lev, double *xi, double *wt, int *nfour, int *nwave) {
  int kl,l,m,n,indd,mmax,jj,maxfxw;
  int nexte=0,nexto=0,next=0,nalpha,nalpha2,j,mm;
  double rb,lscale,x,u,v,w,test1,test2,alpha,halpha,theta;
  double fact[120],cnm[4*mp+1][4*mp+1];
  std::complex<double> ima(0,1);

// YHFRMINI
  fact[0] = 1;
  for( l=1; l<2*mp+2; l++ ) {
    fact[l] = fact[l-1]*l;
  }
  for( l=0; l<mp+1; l++ ) {
    for( m=0; m<=l; m++ ) {
      ytop[l][m] = fact[l-m]/fact[l+m]*(2*l+1);
    }
  }
// YHROTGEN
// BNLCFT
  for( n=0; n<4*mp+1; n++ ) {
    cnm[n][0] = 1;
  }
  for( m=1; m<4*mp+1; m++) {
    cnm[m][m] = 1;
    for( n=m+1; n<4*mp+1; n++ ) {
      cnm[n][m] = cnm[n-1][m]+cnm[n-1][m-1];
    }
  }
  for( m=1; m<4*mp+1; m++) {
    for( n=m+1; n<4*mp+1; n++ ) {
      cnm[n][m] = sqrt(cnm[n][m]);
    }
  }
// FSTRTN
  theta = pi*0.5;
  fstrtn(cnm,rdpi2,theta);
  theta = -pi*0.5;
  fstrtn(cnm,rdmpi2,theta);
  theta = acos(sqrt(3.0)/3);
  fstrtn(cnm,rdsq3,theta);
  theta = acos(-sqrt(3.0)/3);
  fstrtn(cnm,rdmsq3,theta);

// VWGT
  if( mp == 9 ) {
    xi[0]=0.992739967397144e-01;
    xi[1]=0.477256746370494e+00;
    xi[2]=0.105533661382182e+01;
    xi[3]=0.176759343354008e+01;
    xi[4]=0.257342629351470e+01;
    xi[5]=0.344824339201582e+01;
    xi[6]=0.437680983554726e+01;
    xi[7]=0.534895757205460e+01;
    xi[8]=0.635765785313374e+01;
    wt[0]=0.247764418190083e+00;
    wt[1]=0.491885665004643e+00;
    wt[2]=0.653787491376778e+00;
    wt[3]=0.764330384087840e+00;
    wt[4]=0.843761805656281e+00;
    wt[5]=0.904458839850982e+00;
    wt[6]=0.953786131368334e+00;
    wt[7]=0.996702616132185e+00;
    wt[8]=0.104294227302526e+01;
    nfour[0] = 2;
    nfour[1] = 4;
    nfour[2] = 4;
    nfour[3] = 6;
    nfour[4] = 6;
    nfour[5] = 4;
    nfour[6] = 6;
    nfour[7] = 4;
    nfour[8] = 2;
    nwave[0] = 4;
    nwave[1] = 8;
    nwave[2] =12;
    nwave[3] =16;
    nwave[4] =20;
    nwave[5] =20;
    nwave[6] =24;
    nwave[7] = 8;
    nwave[8] = 2;
  } else if ( mp == 18 ) {
    xi[0]=0.527885276611776e-01;
    xi[1]=0.269498598389312e+00;
    xi[2]=0.632203531746893e+00;
    xi[3]=0.111307564277608e+01;
    xi[4]=0.168939496140213e+01;
    xi[5]=0.234376200469530e+01;
    xi[6]=0.306269982907806e+01;
    xi[7]=0.383562941265296e+01;
    xi[8]=0.465424734321562e+01;
    xi[9]=0.551209386593581e+01;
    xi[10]=0.64042126837727e+01;
    xi[11]=0.73268800190617e+01;
    xi[12]=0.82774009925823e+01;
    xi[13]=0.92539718060248e+01;
    xi[14]=0.10255602723746e+02;
    xi[15]=0.11282088297877e+02;
    xi[16]=0.12334067909676e+02;
    xi[17]=0.13414920240172e+02;
    wt[0]=0.134382659143352e+00;
    wt[1]=0.294577527273954e+00;
    wt[2]=0.426078193611486e+00;
    wt[3]=0.531892207765499e+00;
    wt[4]=0.617873062455385e+00;
    wt[5]=0.688631560789050e+00;
    wt[6]=0.747490993814261e+00;
    wt[7]=0.796991927185999e+00;
    wt[8]=0.839174543869975e+00;
    wt[9]=0.875700922837453e+00;
    wt[10]=0.90792943590067e+00;
    wt[11]=0.93698393742461e+00;
    wt[12]=0.96382546688788e+00;
    wt[13]=0.98932985769673e+00;
    wt[14]=0.10143828459791e+01;
    wt[15]=0.10400365437416e+01;
    wt[16]=0.10681548926956e+01;
    wt[17]=0.11090758097553e+01;
    nfour[0]  = 4;
    nfour[1]  = 6;
    nfour[2]  = 6;
    nfour[3]  = 8;
    nfour[4]  = 8;
    nfour[5]  = 8;
    nfour[6]  =10;
    nfour[7]  =10;
    nfour[8]  =10;
    nfour[9]  =10;
    nfour[10] =12;
    nfour[11] =12;
    nfour[12] =12;
    nfour[13] =12;
    nfour[14] =12;
    nfour[15] =12;
    nfour[16] = 8;
    nfour[17] = 2;
    nwave[0]  = 6;
    nwave[1]  = 8;
    nwave[2]  =12;
    nwave[3]  =16;
    nwave[4]  =20;
    nwave[5]  =26;
    nwave[6]  =30;
    nwave[7]  =34;
    nwave[8]  =38;
    nwave[9]  =44;
    nwave[10] =48;
    nwave[11] =52;
    nwave[12] =56;
    nwave[13] =60;
    nwave[14] =60;
    nwave[15] =52;
    nwave[16] = 4;
    nwave[17] = 2;
  } else {
    printf("mp must be 9 or 18 for Helmholtz\n");
    exit(1);
  }
  for( l=0; l<mp; l++ ) {
    test1 = xi[l];
    test2 = sqrt(test1*test1+2*test1*lambda*rd);
    indd = l;
    mmax = nwave[l];
    for( jj=l; jj<mp; jj++ ) {
      if( test2 <= xi[jj] ) {
        indd = jj;
        goto l1001;
      } else {
        mmax = std::max(mmax,nwave[jj]);
      }
    }
l1001: nwave[l] = std::max(mmax,nwave[indd]);
  }

  maxfour = 0;
  maxwave = 0;
  maxfxw = 0;
  nmax = 0;
  for( l=0; l<mp; l++ ) {
    maxfour += nfour[l];
    if( nfour[l] > nmax ) nmax = nfour[l];
    maxwave += nwave[l];
    maxfxw += nfour[l]*nwave[l];
  }
  maxwave = maxwave/2;
  maxexp = std::max(maxfour,maxwave)+1;

// YMKFEXP
  fexpe = new std::complex<double> [maxfxw];
  fexpo = new std::complex<double> [maxfxw];
  fexpb = new std::complex<double> [2*maxfxw];
  for( l=0; l<mp; l++ ) {
    nalpha = nwave[l];
    nalpha2 = nalpha/2;
    halpha = 2*pi/nalpha;
    for( j=0; j<nalpha2; j++ ) {
      alpha = j*halpha;
      for( mm=1; mm<nfour[l]; mm+=2 ) {
        fexpe[nexte] = exp(ima*std::complex<double>(mm*alpha));
        nexte++;
      }
      for( mm=2; mm<nfour[l]; mm+=2 ) {
        fexpo[nexto] = exp(ima*std::complex<double>(mm*alpha));
        nexto++;
      }
    }
    for( mm=2; mm<nfour[l]; mm+=2 ) {
      for( j=0; j<nalpha2; j++ ) {
        alpha = j*halpha;
        fexpb[next] = exp(-ima*std::complex<double>(mm*alpha));
        next++;
      }
    }
    for( mm=1; mm<nfour[l]; mm+=2 ) {
      for( j=0; j<nalpha2; j++ ) {
        alpha = j*halpha;
        fexpb[next] = exp(-ima*std::complex<double>(mm*alpha));
        next++;
      }
    }
  }

  kl = 1 << lev-1;
  rb = rd/kl*scale*0.5;
  lscale = rd/kl*lambda*0.5;

// YRLSCINI
  for( l=0; l<mp; l++ ) {
    x = xi[l]/lscale+1;
    u = sqrt(x*x-1)*rb;
    v = rb*x;
    w = rb*rb;
    m2e[l][0][0] = 1;
    for( m=0; m<mp+1; m++ ) {
      if( m > 0 ) m2e[l][m][m] = m2e[l][m-1][m-1]*u*(2*m-1);
      if( m < mp ) m2e[l][m][m+1] = (2*m+1)*v*m2e[l][m][m];
      for( n=m+2; n<mp+1; n++ ) {
        m2e[l][m][n] = ((2*n-1)*v*m2e[l][m][n-1]-(n+m-1)*w*m2e[l][m][n-2])/(n-m);
      }
    }
  }

// YMKEXPS
  xs = new std::complex<double> [maxwave][3];
  ys = new std::complex<double> [maxwave][3];
  zs = new double [maxwave][3];
  int ntot = 0,ic;
  double w1,w2,hu;
  for( l=0; l<mp; l++ ) {
    w1 = xi[l]+lscale;
    w2 = sqrt( xi[l]*(xi[l]+lscale*2) );
    hu = 2*pi/nwave[l];
    for( m=0; m<nwave[l]/2; m++ ) {
      u = m*hu;
      ic = ntot+m;
      xs[ic][0] = exp(ima*w2*cos(u));
      xs[ic][1] = xs[ic][0]*xs[ic][0];
      xs[ic][2] = xs[ic][1]*xs[ic][0];
      ys[ic][0] = exp(ima*w2*sin(u));
      ys[ic][1] = ys[ic][0]*ys[ic][0];
      ys[ic][2] = ys[ic][1]*ys[ic][0];
      zs[ic][0] = exp(-w1);
      zs[ic][1] = zs[ic][0]*zs[ic][0];
      zs[ic][2] = zs[ic][1]*zs[ic][0];
    }
    ntot += nwave[l]/2;
  }

}
