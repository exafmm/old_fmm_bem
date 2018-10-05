#define MAIN
#include "../include/constants.h"
#include "../include/arrays.h"
#undef MAIN

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi,*vi;
extern float *xj,*yj,*zj,*gxj,*gyj,*gzj,*vj;
extern float *vo;
extern int *nfn,*na,*nb,*nc,*nd;

extern void memoryuse();
extern void memoryfree();
extern void G_dir(int, int, int, int);
extern void Gni_dir(int, int, int, int);
extern void Gnj_dir(int, int, int, int);
extern void alloc(int);
extern void dealloc(int);
#ifdef FMM
extern void nlevel(int);
extern void boxallocate(int, int, int, int);
extern void boxpart(int, int);
extern void sorti(int&);
extern void sortj(int&);
extern void fmm(int, int, int);
extern void unsorti(int&);
extern void unsortj(int&);
#else
extern void tree(int&, int&, int);
#endif

int main(int argc, char *argv[]){
  int neq,ncheck,i,it,np,ni,nj;
  double td,tf,errn,errd,erre;
  std::fstream fid1,fid4;

  neq = 0;
  ncheck = 0;
  umem = 0;
  nprocs = 1;
  myrank = 0;

  xi = new float [npmax];
  yi = new float [npmax];
  zi = new float [npmax];
  gxi = new float [npmax];
  gyi = new float [npmax];
  gzi = new float [npmax];
  vi = new float [npmax];
  xj = new float [npmax];
  yj = new float [npmax];
  zj = new float [npmax];
  gxj = new float [npmax];
  gyj = new float [npmax];
  gzj = new float [npmax];
  vj = new float [npmax];
  vo = new float [npmax];
  nfn = new int [npmax];
  na  = new int [npmax];
  nb  = new int [npmax];
  nc  = new int [npmax];
  nd  = new int [npmax];
  nbi = new int [npmax];
  nbj = new int [npmax];
  sortd = new float [npmax];
  tfmm = new double [100];
  mem = npmax*20*4+100*8;
  memoryuse();

  fid1.open("test.dat",std::ios::out|std::ios::binary);
  fid1.write((char *)(&nprocs),sizeof(int));

  for( i=0; i<npmax; i++ ) {
    xi[i] = rand()/(float) RAND_MAX*2*pi-pi;
    yi[i] = rand()/(float) RAND_MAX*2*pi-pi;
    zi[i] = rand()/(float) RAND_MAX*2*pi-pi;
//    xi[i] = rand()/(float) RAND_MAX*pi;
//    yi[i] = rand()/(float) RAND_MAX*pi;
//    zi[i] = rand()/(float) RAND_MAX*pi;
    gxi[i] = rand()/(float) RAND_MAX;
    gyi[i] = rand()/(float) RAND_MAX;
    gzi[i] = rand()/(float) RAND_MAX;
    xj[i] = rand()/(float) RAND_MAX*2*pi-pi;
    yj[i] = rand()/(float) RAND_MAX*2*pi-pi;
    zj[i] = rand()/(float) RAND_MAX*2*pi-pi;
//    xj[i] = rand()/(float) RAND_MAX*pi-pi;
//    yj[i] = rand()/(float) RAND_MAX*pi-pi;
//    zj[i] = rand()/(float) RAND_MAX*pi-pi;
    gxj[i] = rand()/(float) RAND_MAX;
    gyj[i] = rand()/(float) RAND_MAX;
    gzj[i] = rand()/(float) RAND_MAX;
    vj[i] = rand()/(float) RAND_MAX;
  }

  xmin = -pi;
  ymin = -pi;
  zmin = -pi;
  rd = 2*pi;

  for( it=0; it<25; it++ ) {
    np = int(pow(10,(it+32)/8.0));
    ni = np;
    nj = np;
    printf("N = %d\n",np);

#ifdef FMM
    for( i=0; i<100; i++ ) tfmm[i] = 0;
    nlevel(nj);
    boxpart(ni,nj);
    sorti(ni);
    sortj(nj);
    boxallocate(0,ni,0,nj);
    alloc(1);
    fmm(ni,nj,neq);
    unsorti(ni);
    unsortj(nj);
    dealloc(1);
    tf = 0;
    for ( i=0; i<100; i++ ) tf += tfmm[i];
    printf("fmm    : %g\n",tf);
#else
    tic = get_time();
    tree(ni,nj,neq);
    toc = get_time();
    tf = toc-tic;
    printf("tree   : %g\n",tf);
#endif
    for( i=0; i<ni; i++ ) {
      gxj[i] = gxi[i];
      gyj[i] = gyi[i];
      gzj[i] = gzi[i];
      vo[i] = vi[i];
    }

    alloc(0);
    tic = get_time();
    if( neq == 0 ) {
      G_dir(0,ni-1,0,nj-1);
    } else if( neq == 1 ) {
      Gni_dir(0,ni-1,0,nj-1);
    } else {
      Gnj_dir(0,ni-1,0,nj-1);
    }
    toc = get_time();
    td = toc-tic;
    dealloc(0);
    printf("direct : %g\n",td);

    fid1.write((char *)(&np),sizeof(int));
    for( i=0; i<100; i++ ) fid1.write((char *)(&tfmm[i]),sizeof(double));
    errn = 0;
    for( i=0; i<ni; i++ ) {
      errd = (gxi[i]-gxj[i])*(gxi[i]-gxj[i])+(gyi[i]-gyj[i])*(gyi[i]-gyj[i])+(gzi[i]-gzj[i])*(gzi[i]-gzj[i]);
      erre = gxj[i]*gxj[i]+gyj[i]*gyj[i]+gzj[i]*gzj[i];
//      errd = (vi[i]-vo[i])*(vi[i]-vo[i]);
//      erre = vo[i]*vo[i];
      errn += errd/erre/ni;
    }
    errn = sqrt(errn);
    printf("error  : %g\n\n",errn);
    fid1.write((char *)(&errn),sizeof(double));
  }
  fid1.close();
  delete[] xi;
  delete[] yi;
  delete[] zi;
  delete[] gxi;
  delete[] gyi;
  delete[] gzi;
  delete[] vi;
  delete[] xj;
  delete[] yj;
  delete[] zj;
  delete[] gxj;
  delete[] gyj;
  delete[] gzj;
  delete[] vj;
  delete[] vo;
  delete[] nfn;
  delete[] na;
  delete[] nb;
  delete[] nc;
  delete[] nd;
  delete[] nbi;
  delete[] nbj;
  delete[] sortd;
  delete[] tfmm;
  mem = npmax*20*4+100*8;
  memoryfree();
};
