#include "mpi.h"
#define MAIN
#include "../include/constants.h"
#include "../include/arrays.h"
#undef MAIN

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi,*vi;
extern float *xj,*yj,*zj,*gxj,*gyj,*gzj,*vj;
extern float *gxo,*gyo,*gzo,*vo;
extern int *nfn,*na,*nb,*nc,*nd;
extern double (*tall)[100];

extern void memoryuse();
extern void memoryfree();
extern void G_dir(int, int, int, int);
extern void Gni_dir(int, int, int, int);
extern void Gnj_dir(int, int, int, int);
extern void alloc(int);
extern void dealloc(int);
extern void mpirange(int, int, int&, int&);
extern void mpishift(int, int&, float*);
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
extern void MP_initialize(int*, char***);
extern void MP_end();
extern void tree(int&, int&, int);
#endif

int main(int argc, char **argv){
  int neq,ncheck,mi,mj,mjd,j,i,it,np,ni,nj,ista,iend,jsta,jend;
  double td,tf,errn,errd,erre;
  std::fstream fid1;

  neq = 0;
  ncheck = 0;
  umem = 0;
#ifdef FMM
  MPI_Init(&argc,&argv);
#else
  MP_initialize(&argc,&argv);
#endif
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

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
  gxo = new float [npmax];
  gyo = new float [npmax];
  gzo = new float [npmax];
  vo = new float [npmax];
  vq = new float [npmax];
  nfn = new int [npmax];
  na  = new int [10*npmax];
  nb  = new int [10*npmax];
  nc  = new int [10*npmax];
  nd  = new int [10*npmax];
  nbi = new int [npmax];
  isort = new int [npmax];
  isdsp = new int [nprocs];
  iscnt = new int [nprocs];
  irdsp = new int [nprocs];
  ircnt = new int [nprocs];
  nbj = new int [npmax];
  jsort = new int [npmax];
  jsdsp = new int [nprocs];
  jscnt = new int [nprocs];
  jrdsp = new int [nprocs];
  jrcnt = new int [nprocs];
  sortd = new float [npmax];
  nbuf = new int [npmax];
  nbufd = new int [npmax];
  fbuf = new float [npmax];
  fbufd = new float [npmax];
  dbuf = new double [npmax];
  dbufd = new double [npmax];
  cbuf = new std::complex<float> [npmax];
  cbufd = new std::complex<float> [npmax];
  fsend = new float [npmax];
  frecv = new float [npmax];
  tfmm = new double [100];
  tall = new double [nprocs][100];
  mem = npmax*35*4+npmax*4*8+nprocs*8*8+nprocs*100*8+100*8;
  memoryuse();

  if( myrank == 0 ) {
    fid1.open("test.dat",std::ios::out|std::ios::binary);
    fid1.write((char *)(&nprocs),sizeof(int));
  }

  xmin = -pi;
  ymin = -pi;
  zmin = -pi;
  rd = 2*pi;

  for( it=0; it<25; it++ ) {
    np = int(pow(10,(it+40)/8.0));
    ni = np;
    nj = np;
    if( myrank == 0 ) printf("N = %d\n",np);

    mpirange(0,ni-1,ista,iend);
    mpirange(0,nj-1,jsta,jend);
    mi = iend-ista+1;
    mj = jend-jsta+1;

    srand(myrank);
    for( i=0; i<mi; i++ ) {
      xi[i] = rand()/(float) RAND_MAX*2*pi-pi;
      yi[i] = rand()/(float) RAND_MAX*2*pi-pi;
      zi[i] = rand()/(float) RAND_MAX*2*pi-pi;
//      xi[i] = rand()/(float) RAND_MAX*pi;
//      yi[i] = rand()/(float) RAND_MAX*pi;
//      zi[i] = rand()/(float) RAND_MAX*pi;
      gxi[i] = rand()/(float) RAND_MAX;
      gyi[i] = rand()/(float) RAND_MAX;
      gzi[i] = rand()/(float) RAND_MAX;
    }
    for( i=0; i<mj; i++ ) {
      xj[i] = rand()/(float) RAND_MAX*2*pi-pi;
      yj[i] = rand()/(float) RAND_MAX*2*pi-pi;
      zj[i] = rand()/(float) RAND_MAX*2*pi-pi;
//      xj[i] = rand()/(float) RAND_MAX*pi-pi;
//      yj[i] = rand()/(float) RAND_MAX*pi-pi;
//      zj[i] = rand()/(float) RAND_MAX*pi-pi;
      gxj[i] = rand()/(float) RAND_MAX;
      gyj[i] = rand()/(float) RAND_MAX;
      gzj[i] = rand()/(float) RAND_MAX;
      vj[i] = rand()/(float) RAND_MAX;
    }

#ifdef FMM
    for( i=0; i<100; i++ ) tfmm[i] = 0;
    nlevel(mj);
    boxpart(mi,mj);
    sorti(mi);
    sortj(mj);
    boxallocate(0,mi,0,mj);
    alloc(1);
    fmm(mi,mj,neq);
    unsorti(mi);
    unsortj(mj);
    dealloc(1);
    tf = 0;
//    for ( i=0; i<100; i++ ) tf += tfmm[i];
    for ( i=10; i<23; i++ ) tf += tfmm[i];
    MPI_Gather(tfmm,100,MPI_DOUBLE,tall,100,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if( myrank == 0 ) printf("fmm    : %g\n",tf);
#else
    MPI_Barrier(MPI_COMM_WORLD);
    tic = MPI_Wtime();
    tree(mi,mj,neq);
    MPI_Barrier(MPI_COMM_WORLD);
    toc = MPI_Wtime();
    tf = toc-tic;
    if( myrank == 0 ) printf("tree   : %g\n",tf);
#endif
    for( i=0; i<mi; i++ ) {
      gxj[i] = gxi[i];
      gyj[i] = gyi[i];
      gzj[i] = gzi[i];
      vq[i] = vi[i];
    }

#if 0
    alloc(0);
    MPI_Barrier(MPI_COMM_WORLD);
    tic = MPI_Wtime();
    for( i=0; i<mi; i++ ) {
      gxo[i] = 0;
      gyo[i] = 0;
      gzo[i] = 0;
      vo[i] = 0;
    }
    for( j=0; j<nprocs; j++ ) {
      mpishift(mj,mjd,xj);
      mpishift(mj,mjd,yj);
      mpishift(mj,mjd,zj);
      mpishift(mj,mjd,gxj);
      mpishift(mj,mjd,gyj);
      mpishift(mj,mjd,gzj);
      mpishift(mj,mjd,vj);
      mj = mjd;
      if( neq == 0 ) {
        G_dir(0,mi-1,0,mj-1);
      } else if( neq == 1 ) {
        Gni_dir(0,mi-1,0,mj-1);
      } else  {
        Gnj_dir(0,mi-1,0,mj-1);
      }
      for( i=0; i<mi; i++ ) {
        gxo[i] += gxi[i];
        gyo[i] += gyi[i];
        gzo[i] += gzi[i];
        vo[i] += vi[i];
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    toc = MPI_Wtime();
    td = toc-tic;
    dealloc(0);
#endif
    if( myrank == 0 ) printf("direct : %g\n",td);

    errn = 0;
    for( i=0; i<mi; i++ ) {
      errd = (gxo[i]-gxj[i])*(gxo[i]-gxj[i])+(gyo[i]-gyj[i])*(gyo[i]-gyj[i])+(gzo[i]-gzj[i])*(gzo[i]-gzj[i]);
      erre = gxj[i]*gxj[i]+gyj[i]*gyj[i]+gzj[i]*gzj[i];
//      errd = (vq[i]-vo[i])*(vq[i]-vo[i]);
//      erre = vo[i]*vo[i];
      errn += errd/erre/ni;
    }
    MPI_Reduce(&errn,&errd,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    if( myrank == 0 ) {
      fid1.write((char *)(&np),sizeof(int));
      for( j=0; j<nprocs; j++ ) {
        for( i=0; i<100; i++ ) {
          fid1.write((char *)(&tall[j][i]),sizeof(double));
        }
      }
      errn = sqrt(errd);
      printf("error  : %g\n\n",errn);
      fid1.write((char *)(&errn),sizeof(double));
    }
  }

  if( myrank == 0 ) {
    fid1.close();
  }
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
  delete[] gxo;
  delete[] gyo;
  delete[] gzo;
  delete[] vo;
  delete[] vq;
  delete[] nfn;
  delete[] na;
  delete[] nb;
  delete[] nc;
  delete[] nd;
  delete[] nbi;
  delete[] isort;
  delete[] isdsp;
  delete[] iscnt;
  delete[] irdsp;
  delete[] ircnt;
  delete[] nbj;
  delete[] jsort;
  delete[] jsdsp;
  delete[] jscnt;
  delete[] jrdsp;
  delete[] jrcnt;
  delete[] sortd;
  delete[] nbuf;
  delete[] nbufd;
  delete[] fbuf;
  delete[] fbufd;
  delete[] dbuf;
  delete[] dbufd;
  delete[] cbuf;
  delete[] cbufd;
  delete[] fsend;
  delete[] frecv;
  delete[] tfmm;
  delete[] tall;
  mem = npmax*35*4+nprocs*8*8+nprocs*100*8+100*8;
  memoryfree();
#ifdef FMM
  MPI_Finalize();
#else
  MP_end();
#endif
};
