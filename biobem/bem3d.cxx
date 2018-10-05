#include "../include/constants.h"

extern float (*ptr)[3];
extern int (*ntr)[6],(*nne)[20];
extern float (*vna)[3],(*vta)[3],(*vsa)[3];
extern float (*p1)[13],(*p2)[13],(*p3)[13],(*p4)[13],(*p5)[13],(*p6)[13];
extern float (*hsg)[13],*are,*crv,*xiq,*etq,*wq;

void gauss_trgl() {
  double al,be,ga,de,o1,o2,qa,rh,ru,o3,o4;

  if( ngl == 1 ) {
    xiq[0] = 1.0/3;
    etq[0] = 1.0/3;
     wq[0] = 1.0;
  } else if( ngl == 3 ) {
    xiq[0] = 1.0/6;
    etq[0] = 1.0/6;
     wq[0] = 1.0/3;
    
    xiq[1] = 2.0/3;
    etq[1] = 1.0/6;
     wq[1] = wq[0];

    xiq[2] = 1.0/6;
    etq[2] = 2.0/3;
     wq[2] = wq[0];
  } else if( ngl == 4 ) {
    xiq[0] = 1.0/3;
    etq[0] = 1.0/3;
     wq[0] = -27.0/48;

    xiq[1] = 1.0/5;
    etq[1] = 1.0/5;
     wq[1] = 25.0/48;

    xiq[2] = 3.0/5;
    etq[2] = 1.0/5;
     wq[2] = 25.0/48;

    xiq[3] = 1.0/5;
    etq[3] = 3.0/5;
     wq[3] = 25.0/48;
  } else if( ngl == 6 ) {
    al = 0.816847572980459;
    be = 0.445948490915965;
    ga = 0.108103018168070;
    de = 0.091576213509771;
    o1 = 0.109951743655322;
    o2 = 0.223381589678011;

    xiq[0] = de;
    xiq[1] = al;
    xiq[2] = de;
    xiq[3] = be;
    xiq[4] = ga;
    xiq[5] = be;

    etq[0] = de;
    etq[1] = de;
    etq[2] = al;
    etq[3] = be;
    etq[4] = be;
    etq[5] = ga;

    wq[0] = o1;
    wq[1] = o1;
    wq[2] = o1;
    wq[3] = o2;
    wq[4] = o2;
    wq[5] = o2;
  } else if( ngl == 7 ) {
    al = 0.797426958353087;
    be = 0.470142064105115;
    ga = 0.059715871789770;
    de = 0.101286507323456;
    o1 = 0.125939180544827;
    o2 = 0.132394152788506;

    xiq[0] = de;
    xiq[1] = al;
    xiq[2] = de;
    xiq[3] = be;
    xiq[4] = ga;
    xiq[5] = be;
    xiq[6] = 1.0/3;

    etq[0] = de;
    etq[1] = de;
    etq[2] = al;
    etq[3] = be;
    etq[4] = be;
    etq[5] = ga;
    etq[6] = 1.0/3;

    wq[0] = o1;
    wq[1] = o1;
    wq[2] = o1;
    wq[3] = o2;
    wq[4] = o2;
    wq[5] = o2;
    wq[6] = 0.225;
  } else if( ngl == 9 ) {
    al = 0.124949503233232;
    qa = 0.165409927389841;
    rh = 0.797112651860071;
    de = 0.437525248383384;
    ru = 0.037477420750088;
    o1 = 0.205950504760887;
    o2 = 0.063691414286223;

    xiq[0] = de;
    xiq[1] = al;
    xiq[2] = de;
    xiq[3] = qa;
    xiq[4] = ru;
    xiq[5] = rh;
    xiq[6] = qa;
    xiq[7] = ru;
    xiq[8] = rh;

    etq[0] = de;
    etq[1] = de;
    etq[2] = al;
    etq[3] = ru;
    etq[4] = qa;
    etq[5] = qa;
    etq[6] = rh;
    etq[7] = rh;
    etq[8] = ru;

    wq[0] = o1;
    wq[1] = o1;
    wq[2] = o1;
    wq[3] = o2;
    wq[4] = o2;
    wq[5] = o2;
    wq[6] = o2;
    wq[7] = o2;
    wq[8] = o2;
  } else if( ngl == 12 ) {
    al = 0.873821971016996;
    be = 0.249286745170910;
    ga = 0.501426509658179;
    de = 0.063089014491502;
    rh = 0.636502499121399;
    qa = 0.310352451033785;
    ru = 0.053145049844816;
    o1 = 0.050844906370207;
    o2 = 0.116786275726379;
    o3 = 0.082851075618374;
 
    xiq[0]  = de;
    xiq[1]  = al;
    xiq[2]  = de;
    xiq[3]  = be;
    xiq[4]  = ga;
    xiq[5]  = be;
    xiq[6]  = qa;
    xiq[7]  = ru;
    xiq[8]  = rh;
    xiq[9]  = qa;
    xiq[10] = ru;
    xiq[11] = rh;

    etq[0]  = de;
    etq[1]  = de;
    etq[2]  = al;
    etq[3]  = be;
    etq[4]  = be;
    etq[5]  = ga;
    etq[6]  = ru;
    etq[7]  = qa;
    etq[8]  = qa;
    etq[9]  = rh;
    etq[10] = rh;
    etq[11] = ru;

    wq[0]  = o1;
    wq[1]  = o1;
    wq[2]  = o1;
    wq[3]  = o2;
    wq[4]  = o2;
    wq[5]  = o2;
    wq[6]  = o3;
    wq[7]  = o3;
    wq[8]  = o3;
    wq[9]  = o3;
    wq[10] = o3;
    wq[11] = o3;
  } else if( ngl == 13 ) {
    al = 0.479308067841923;
    be = 0.065130102902216;
    ga = 0.869739794195568;
    de = 0.260345966079038;
    rh = 0.638444188569809;
    qa = 0.312865496004875;
    ru = 0.048690315425316;
    o1 = 0.175615257433204;
    o2 = 0.053347235608839;
    o3 = 0.077113760890257;
    o4 =-0.149570044467670;

    xiq[0]  = de;
    xiq[1]  = al;
    xiq[2]  = de;
    xiq[3]  = be;
    xiq[4]  = ga;
    xiq[5]  = be;
    xiq[6]  = qa;
    xiq[7]  = ru;
    xiq[8]  = rh;
    xiq[9]  = qa;
    xiq[10] = ru;
    xiq[11] = rh;
    xiq[12] = 1.0/3;

    etq[0]  = de;
    etq[1]  = de;
    etq[2]  = al;
    etq[3]  = be;
    etq[4]  = be;
    etq[5]  = ga;
    etq[6]  = ru;
    etq[7]  = qa;
    etq[8]  = qa;
    etq[9]  = rh;
    etq[10] = rh;
    etq[11] = ru;
    etq[12] = 1.0/3;

    wq[0]  = o1;
    wq[1]  = o1;
    wq[2]  = o1;
    wq[3]  = o2;
    wq[4]  = o2;
    wq[5]  = o2;
    wq[6]  = o3;
    wq[7]  = o3;
    wq[8]  = o3;
    wq[9]  = o3;
    wq[10] = o3;
    wq[11] = o3;
    wq[12] = o4;
  }
}

void interp_p(double x[], double y[], double z[],
              double al, double be, double ga, double xi, double eta, 
              double& ph1, double& ph2, double& ph3, double& ph4, double& ph5, double& ph6,
              double& dxdxi, double& dydxi, double& dzdxi,
              double& dxdet, double& dydet, double& dzdet,
              double& vnx, double& vny, double& vnz, double& hxi, double& het, double& hsd) {
  double alc,bec,gac,alalc,bebec,gagac,dph1,dph2,dph3,dph4,dph5,dph6,pph1,pph2,pph3,pph4,pph5,pph6;

  alc = 1.0-al;
  bec = 1.0-be;
  gac = 1.0-ga;
  alalc = al*alc;
  bebec = be*bec;
  gagac = ga*gac;

// Shape functions

  ph2 = xi *(xi-al+eta*(al-ga)/gac)/alc;
  ph3 = eta*(eta-be+xi*(be+ga-1.0)/ga)/bec;
  ph4 = xi *(1.0-xi-eta)/alalc;
  ph5 = xi*eta/gagac;
  ph6 = eta*(1.0-xi-eta)/bebec;
  ph1 = 1.0-ph2-ph3-ph4-ph5-ph6;

// Xi derivatives of phi

  dph2 = (2.0*xi-al+eta*(al-ga)/gac)/alc;
  dph3 = eta*(be+ga-1.0)/(ga*bec);
  dph4 = (1.0-2.0*xi-eta)/alalc;
  dph5 = eta/gagac;
  dph6 = -eta/bebec;
  dph1 = -dph2-dph3-dph4-dph5-dph6;

// Compute dx/dxi from xi derivatives of phi

  dxdxi = x[0]*dph1+x[1]*dph2+x[2]*dph3+x[3]*dph4+x[4]*dph5+x[5]*dph6;
  dydxi = y[0]*dph1+y[1]*dph2+y[2]*dph3+y[3]*dph4+y[4]*dph5+y[5]*dph6;
  dzdxi = z[0]*dph1+z[1]*dph2+z[2]*dph3+z[3]*dph4+z[4]*dph5+z[5]*dph6;

// Eta derivatives of phi

  pph2 = xi*(al-ga)/(alc*gac);
  pph3 = (2.0*eta-be+xi*(be+ga-1.0)/ga)/bec;
  pph4 = -xi/alalc;
  pph5 = xi/gagac;
  pph6 = (1.0-xi-2.0*eta)/bebec;
  pph1 = -pph2-pph3-pph4-pph5-pph6;

// Compute dx/deta from eta derivatives of phi

  dxdet = x[0]*pph1+x[1]*pph2+x[2]*pph3+x[3]*pph4+x[4]*pph5+x[5]*pph6;
  dydet = y[0]*pph1+y[1]*pph2+y[2]*pph3+y[3]*pph4+y[4]*pph5+y[5]*pph6;
  dzdet = z[0]*pph1+z[1]*pph2+z[2]*pph3+z[3]*pph4+z[4]*pph5+z[5]*pph6;

// Normal vector    vn = (dxdxi)x(dxdeta)
// Surface metric   hsd = norm(vn)

  vnx = dydxi*dzdet-dydet*dzdxi;
  vny = dzdxi*dxdet-dzdet*dxdxi;
  vnz = dxdxi*dydet-dxdet*dydxi;
  hsd  = sqrt(vnx*vnx+vny*vny+vnz*vnz);

// Normalizations

//  hxi = sqrt(dxdxi*dxdxi+dydxi*dydxi+dzdxi*dzdxi);
//  dxdxi = dxdxi/hxi;
//  dydxi = dydxi/hxi;
//  dzdxi = dzdxi/hxi;
//  het = sqrt(dxdet*dxdet+dydet*dydet+dzdet*dzdet)
//  dxdet = dxdet/het;
//  dydet = dydet/het;
//  dzdet = dzdet/het;
  vnx = vnx/hsd;
  vny = vny/hsd;
  vnz = vnz/hsd;
}

void elm_geom() {
  int i,j,m;
  double d42,d41,d63,d61,d52,d53,al,be,ga,alc,bec,gac,xi,eta;
  double ph1,ph2,ph3,ph4,ph5,ph6,hxi,het,hsd,cf;
  double bvx1,bvy1,bvz1,bvx2,bvy2,bvz2,bvx3,bvy3,bvz3,crvx,crvy,crvz;
  double dxdxi,dydxi,dzdxi,dxdei,dydei,dzdei,vnx,vny,vnz,par,sabs;
  double xxi[6],eet[6],dxdx[6],dydx[6],dzdx[6],dxde[6],dyde[6],dzde[6];
  double vx[6],vy[6],vz[6],x[6],y[6],z[6];

  for( i=0; i<nwn; i++ ) {
    vna[i][0] = 0;
    vna[i][1] = 0;
    vna[i][2] = 0;

    vta[i][0] = 0;
    vta[i][1] = 0;
    vta[i][2] = 0;

    vsa[i][0] = 0;
    vsa[i][1] = 0;
    vsa[i][2] = 0;
  }
  for( i=0; i<nwe; i++ ) {
    x[0] = ptr[ntr[i][0]][0];
    x[1] = ptr[ntr[i][1]][0];
    x[2] = ptr[ntr[i][2]][0];
    x[3] = ptr[ntr[i][3]][0];
    x[4] = ptr[ntr[i][4]][0];
    x[5] = ptr[ntr[i][5]][0];

    y[0] = ptr[ntr[i][0]][1];
    y[1] = ptr[ntr[i][1]][1];
    y[2] = ptr[ntr[i][2]][1];
    y[3] = ptr[ntr[i][3]][1];
    y[4] = ptr[ntr[i][4]][1];
    y[5] = ptr[ntr[i][5]][1];

    z[0] = ptr[ntr[i][0]][2];
    z[1] = ptr[ntr[i][1]][2];
    z[2] = ptr[ntr[i][2]][2];
    z[3] = ptr[ntr[i][3]][2];
    z[4] = ptr[ntr[i][4]][2];
    z[5] = ptr[ntr[i][5]][2];

    d42 = sqrt((x[3]-x[1])*(x[3]-x[1])+(y[3]-y[1])*(y[3]-y[1])+(z[3]-z[1])*(z[3]-z[1]));
    d41 = sqrt((x[3]-x[0])*(x[3]-x[0])+(y[3]-y[0])*(y[3]-y[0])+(z[3]-z[0])*(z[3]-z[0]));
    d63 = sqrt((x[5]-x[2])*(x[5]-x[2])+(y[5]-y[2])*(y[5]-y[2])+(z[5]-z[2])*(z[5]-z[2]));
    d61 = sqrt((x[5]-x[0])*(x[5]-x[0])+(y[5]-y[0])*(y[5]-y[0])+(z[5]-z[0])*(z[5]-z[0]));
    d52 = sqrt((x[4]-x[1])*(x[4]-x[1])+(y[4]-y[1])*(y[4]-y[1])+(z[4]-z[1])*(z[4]-z[1]));
    d53 = sqrt((x[4]-x[2])*(x[4]-x[2])+(y[4]-y[2])*(y[4]-y[2])+(z[4]-z[2])*(z[4]-z[2]));

    al = 1.0/(1.0+d42/d41);
    be = 1.0/(1.0+d63/d61);
    ga = 1.0/(1.0+d52/d53);
    alc = 1.0-al;
    bec = 1.0-be;
    gac = 1.0-ga;
    are[i] = 0;

    xxi[0] = 0;
    eet[0] = 0;
    xxi[1] = 1;
    eet[1] = 0;
    xxi[2] = 0;
    eet[2] = 1;
    xxi[3] = al;
    eet[3] = 0;
    xxi[4] = ga;
    eet[4] = gac;
    xxi[5] = 0;
    eet[5] = be;

// Compute shape function & surface metric

    for( j=0; j<ngl; j++ ) {
      xi = xiq[j];
      eta = etq[j];

      interp_p(x,y,z,al,be,ga,xi,eta,ph1,ph2,ph3,ph4,ph5,ph6,
               dxdx[j],dydx[j],dzdx[j],dxde[j],dyde[j],dzde[j],
               vx[j],vy[j],vz[j],hxi,het,hsd);

      cf = hsd*wq[j];
      are[i] += cf;
      hsg[i][j] = cf;

      p1[i][j] = ph1;
      p2[i][j] = ph2;
      p3[i][j] = ph3;
      p4[i][j] = ph4;
      p5[i][j] = ph5;
      p6[i][j] = ph6;
    }
    are[i] *= 0.5;

// Compute normal vector

    for( j=0; j<6; j++ ) {
      xi = xxi[j];
      eta = eet[j];

      interp_p(x,y,z,al,be,ga,xi,eta,ph1,ph2,ph3,ph4,ph5,ph6,
               dxdx[j],dydx[j],dzdx[j],dxde[j],dyde[j],dzde[j],
               vx[j],vy[j],vz[j],hxi,het,hsd);

      m = ntr[i][j];
      vna[m][0] += vx[j];
      vna[m][1] += vy[j];
      vna[m][2] += vz[j];
    }
    if( 0 ) { // Skip mean curvature computation

// Compute mean curvature

// Segment 1-4-2

      bvx1 = vy[0]*dzdx[0]-vz[0]*dydx[0];
      bvy1 = vz[0]*dxdx[0]-vx[0]*dzdx[0];
      bvz1 = vx[0]*dydx[0]-vy[0]*dxdx[0];
      bvx2 = vy[3]*dzdx[3]-vz[3]*dydx[3];
      bvy2 = vz[3]*dxdx[3]-vx[3]*dzdx[3];
      bvz2 = vx[3]*dydx[3]-vy[3]*dxdx[3];
      bvx3 = vy[1]*dzdx[1]-vz[1]*dydx[1];
      bvy3 = vz[1]*dxdx[1]-vx[1]*dzdx[1];
      bvz3 = vx[1]*dydx[1]-vy[1]*dxdx[1];

      crvx = al*bvx1+bvx2+alc*bvx3;
      crvy = al*bvy1+bvy2+alc*bvy3;
      crvz = al*bvz1+bvz2+alc*bvz3;

// Segment 2-5-3

      bvx1 = vy[1]*dzdx[1]-vz[1]*dydx[1];
      bvy1 = vz[1]*dxdx[1]-vx[1]*dzdx[1];
      bvz1 = vx[1]*dydx[1]-vy[1]*dxdx[1];
      bvx2 = vy[4]*dzdx[4]-vz[4]*dydx[4];
      bvy2 = vz[4]*dxdx[4]-vx[4]*dzdx[4];
      bvz2 = vx[4]*dydx[4]-vy[4]*dxdx[4];
      bvx3 = vy[2]*dzdx[2]-vz[2]*dydx[2];
      bvy3 = vz[2]*dxdx[2]-vx[2]*dzdx[2];
      bvz3 = vx[2]*dydx[2]-vy[2]*dxdx[2];

      crvx = crvx-gac*bvx1-bvx2-ga*bvx3;
      crvy = crvy-gac*bvy1-bvy2-ga*bvy3;
      crvz = crvz-gac*bvz1-bvz2-ga*bvz3;

      bvx1 = vy[1]*dzde[1]-vz[1]*dyde[1];
      bvy1 = vz[1]*dxde[1]-vx[1]*dzde[1];
      bvz1 = vx[1]*dyde[1]-vy[1]*dxde[1];
      bvx2 = vy[4]*dzde[4]-vz[4]*dyde[4];
      bvy2 = vz[4]*dxde[4]-vx[4]*dzde[4];
      bvz2 = vx[4]*dyde[4]-vy[4]*dxde[4];
      bvx3 = vy[2]*dzde[2]-vz[2]*dyde[2];
      bvy3 = vz[2]*dxde[2]-vx[2]*dzde[2];
      bvz3 = vx[2]*dyde[2]-vy[2]*dxde[2];

      crvx = crvx+gac*bvx1+bvx2+ga*bvx3;
      crvy = crvy+gac*bvy1+bvy2+ga*bvy3;
      crvz = crvz+gac*bvz1+bvz2+ga*bvz3;

// Segment 3-6-1

      bvx1 = vy[0]*dzde[0]-vz[0]*dyde[0];
      bvy1 = vz[0]*dxde[0]-vx[0]*dzde[0];
      bvz1 = vx[0]*dyde[0]-vy[0]*dxde[0];
      bvx2 = vy[5]*dzde[5]-vz[5]*dyde[5];
      bvy2 = vz[5]*dxde[5]-vx[5]*dzde[5];
      bvz2 = vx[5]*dyde[5]-vy[5]*dxde[5];
      bvx3 = vy[2]*dzde[2]-vz[2]*dyde[2];
      bvy3 = vz[2]*dxde[2]-vx[2]*dzde[2];
      bvz3 = vx[2]*dyde[2]-vy[2]*dxde[2];

      crvx = crvx-be*bvx1-bvx2-bec*bvx3;
      crvy = crvy-be*bvy1-bvy2-bec*bvy3;
      crvz = crvz-be*bvz1-bvz2-bec*bvz3;
      cf = 0.25/are[i];

// Average normal vector at centroid

      xi = 1.0/3;
      eta = 1.0/3;

      interp_p(x,y,z,al,be,ga,xi,eta,ph1,ph2,ph3,ph4,ph5,ph6,
               dxdxi,dydxi,dzdxi,dxdei,dydei,dzdei,
               vnx,vny,vnz,hxi,het,hsd);

      crv[i] = cf*(crvx*vnx+crvy*vny+crvz*vnz);
    }
  }

// Normalize normal vector & compute tangential & binormal vectors

  for( i=0; i<nwn; i++ ) {
    par = (double) nne[i][0];
    vna[i][0] /= par;
    vna[i][1] /= par;
    vna[i][2] /= par;

    par = sqrt(vna[i][0]*vna[i][0]+vna[i][1]*vna[i][1]+vna[i][2]*vna[i][2]);
    vna[i][0] /= par;
    vna[i][1] /= par;
    vna[i][2] /= par;

    sabs = sqrt(vna[i][1]*vna[i][1]+vna[i][2]*vna[i][2]);
    if( sabs > 1e-3 ) {
      vsa[i][0] = 0;
      vsa[i][1] = vna[i][2]/sabs;
      vsa[i][2] = -vna[i][1]/sabs;
    } else {
      vsa[i][0] = 0;
      vsa[i][1] = 1;
      vsa[i][2] = 0;
    }
    vta[i][0] = vsa[i][1]*vna[i][2]-vsa[i][2]*vna[i][1];
    vta[i][1] = vsa[i][2]*vna[i][0]-vsa[i][0]*vna[i][2];
    vta[i][2] = vsa[i][0]*vna[i][1]-vsa[i][1]*vna[i][0];
  }
}
