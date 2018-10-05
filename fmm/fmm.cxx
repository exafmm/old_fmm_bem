#include "../include/constants.h"

extern int *nlbi,*nej,*nfj,*nlbj;
extern float *gxi,*gyi,*gzi,*vi;
extern double *tfmm;

extern void precalc(int&, int);
extern void boxdatai(int, int&, double&);
extern void boxdataj(int, int&, double&);
extern void boxm2mi(int&, int);
extern void boxm2mj(int&, int);
extern void ifbox(int, int, int, int, int, int&, int, int);
extern void G_p2p(int, int);
extern void Gni_p2p(int, int);
extern void Gnj_p2p(int, int);
extern void Helmholtz_G_p2p(int, int);
extern void G_p2m(int, int, int);
extern void Gn_p2m(int, int, int);
extern void Helmholtz_G_p2m(int, int, int);
extern void m2m(int, int, int, int, int);
extern void Helmholtz_m2m(int, int, int, int, int);
extern void m2l(int, int, int, int, int, int);
extern void Helmholtz_m2l(int, int, int, int, int, int);
extern void l2l(int, int, int, int);
extern void G_l2p(int, int, int);
extern void Gn_l2p(int, int, int);
extern void Helmholtz_G_l2p(int, int, int);

void fmm(int mi, int mj, int neq){
  int lev,lmin,i,nmp,lbi,lbj,lbjr,lbjo;
  double rb;

  lev = lmax;
  lmin = 2;
  nlbi[lev-1] = 0;
  nlbj[lev-1] = 0;

  precalc(nmp,mp);
  boxdatai(mi,lbi,rb);
  boxdataj(mj,lbj,rb);

// Step 1. P2P-summation @ lev lmax

  ifbox(mi,mj,nmp,lbi,lbj,lbjr,lev,0);

  for( i=0; i<mi; i++ ) vi[i] = 0;
  if( neq == 0 ) {
    for( i=0; i<mi; i++ ) {
      gxi[i] = 0;
      gyi[i] = 0;
      gzi[i] = 0;
    }
    G_p2p(lbi,lbj);
  } else if( neq == 1 ) {
    Gni_p2p(lbi,lbj);
  } else if( neq == 2 ) {
    Gnj_p2p(lbi,lbj);
  } else if( neq == 10 ) {
    for( i=0; i<mi; i++ ) {
      gxi[i] = 0;
      gyi[i] = 0;
      gzi[i] = 0;
    }
    Helmholtz_G_p2p(lbi,lbj);
  }

  lev = lmax;

// Step 2. P2M-expansion @ lev lmax

  if( neq == 0 || neq == 1 ) {
    G_p2m(nmp,mp,lbj);
  } else if( neq == 2 ) {
    Gn_p2m(nmp,mp,lbj);
  } else if( neq == 10 ) {
    Helmholtz_G_p2m(nmp,mp,lbj);
  }

// Step 3. M2M-translation @ lev lmax-1 to 2

  if(lmax > lmin) {
    for( lev=lmax-1; lev>=lmin; lev-- ) {
      lbjo = lbj;
      boxm2mi(lbi,lev);
      boxm2mj(lbj,lev);
      if( neq < 10 ) {
        m2m(nmp,mp,lev,lbj,lbjo);
      } else if ( neq < 20 ) {
        Helmholtz_m2m(nmp,mp,lev,lbj,lbjo);
      }
    }
    lev = lmin;
  } else {
    int jj,jb;
    for( i=0; i<nbmax; i++ ) nej[i] = -1;
    for( jj=0; jj<lbj; jj++ ) {
      jb = jj+nlbj[lev-1];
      nej[nfj[jb]] = jj; 
    }
  }

// Step 4. M2L-translation @ lev 2

  ifbox(mi,mj,nmp,lbi,lbj,lbjr,lev,3);

  if( neq < 10 ) {
    m2l(nmp,mp,lbi,lbjr,lev,lbi);
  } else if( neq < 20 ) {
    Helmholtz_m2l(nmp,mp,lbi,lbjr,lev,lbi);
  }

// Step 5. L2L-translation @ lev 3 to lmax

  if( lmax > lmin ) {
    for( lev=lmin+1; lev<=lmax; lev++ ) {
      lbi = nlbi[lev-2]-nlbi[lev-1];
      lbj = nlbj[lev-2]-nlbj[lev-1];

      l2l(nmp,mp,lev,lbi);

      int jj,jb;
      for( i=0; i<nbmax; i++ ) nej[i] = -1;
      for( jj=0; jj<lbj; jj++ ) {
        jb = jj+nlbj[lev-1];
        nej[nfj[jb]] = jj;
      }

// Step 6. M2L-translation @ lev 3 to lmax

      ifbox(mi,mj,nmp,lbi,lbj,lbjr,lev,2);

      if( neq < 10 ) {
        m2l(nmp,mp,lbi,lbjr,lev,0);
      } else if( neq < 20 ) {
        Helmholtz_m2l(nmp,mp,lbi,lbjr,lev,0);
      }

    }
    lev = lmax;
  }

// Step 7. L2P-expansion @ lev lmax

  if( neq == 0 || neq == 2 ) {
    G_l2p(nmp,mp,lbi);
  } else if( neq == 1 ) {
    Gn_l2p(nmp,mp,lbi);
  } else if( neq < 20 ) {
    Helmholtz_G_l2p(nmp,mp,lbi);
  }

}
