#define MAIN
#include "petscsnes.h"
#include "parameters.h"
#include "../include/constants.h"
#include "../include/arrays.h"
#undef MAIN
#define DIR 0

extern void memoryuse();
extern void memoryfree();
extern void mpirange(int, int, int&, int&);
extern void nlevel(int);
extern void boxallocate(int, int, int, int);
extern void alloc(int);
extern void boxpart(int, int);
extern void sorti(int&);
extern void sortj(int&);
extern void unsorti(int&);
extern void unsortj(int&);
extern void dealloc(int);

extern void mpishift(int, int&, float*);
extern void G_dir(int, int, int, int);
extern void Gni_dir(int, int, int, int);
extern void Gnj_dir(int, int, int, int);
extern void fmm(int, int, int);

PetscErrorCode mymatmult(Mat A,Vec X,Vec Y)
{
  PetscInt mi,mj,mjd,ista,iend;
  PetscScalar tic,toc,*x,*y;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = VecSet(Y,0.0);CHKERRQ(ierr);
  ierr = VecGetArray(X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(Y,&y);CHKERRQ(ierr);
  ierr = VecGetOwnershipRange(X,&ista,&iend);CHKERRQ(ierr);
  mi = iend-ista;
  mj = iend-ista;
  for(PetscInt m = 0; m < mj; ++m) vj[m] = x[m]*are[m]/e1;

  tic = get_time();
//  if(myrank==0) std::cout << "kernel start, mi : " << mi << std::endl;
#if DIR
  for(PetscInt l = 0; l < mi; ++l) vo[l] = 0;
  for(PetscInt m = 0; m < nprocs; ++m) {
    mpishift(mj,mjd,xj);
    mpishift(mj,mjd,yj);
    mpishift(mj,mjd,zj);
    mpishift(mj,mjd,vj);
    mj = mjd;
    Gni_dir(0,mi-1,0,mj-1);
    for(PetscInt l = 0; l < mi; ++l) vo[l] += vi[l];
  }
  for(PetscInt l = 0; l < mi; ++l) vi[l] = vo[l];
#else
  fmm(mi,mj,1);
#endif
  toc = get_time();
//  if(myrank==0) std::cout << "kernel end, t : " << toc-tic << std::endl;
/*
  if(myrank==0) {
    for(PetscInt l = 0; l < 28; ++l) {
      std::cout << "tfmm[" << l << "] : " << tfmm[l] << std::endl;
    }
    toc = tfmm[12];
    for(PetscInt l = 17; l < 23; ++l) toc += tfmm[l];
    std::cout << "real kernel, t : " << toc << std::endl;
  }
*/

  for(PetscInt l = 0; l < mi; ++l) y[l] = are[l]*(vi[l]+(e1+e2)/(2*e1*(e1-e2))*x[l]);

  ierr = VecRestoreArray(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(Y,&y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode mypcshell(PC pc,Vec X,Vec Y)
{
  PetscInt mi,ista,iend;
  PetscScalar *x,*y;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = VecSet(Y,0.0);CHKERRQ(ierr);
  ierr = VecGetArray(X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(Y,&y);CHKERRQ(ierr);
  ierr = VecGetOwnershipRange(X,&ista,&iend);CHKERRQ(ierr);
  mi = iend-ista;

  for(PetscInt l = 0; l < mi; ++l) y[l] = x[l]/are[l]/(e1+e2)*(2*e1*(e1-e2));
//  for(PetscInt l = 0; l < mi; ++l) y[l] = x[l];

  ierr = VecRestoreArray(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(Y,&y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

int main(int argc, char *argv[]) {
  Mat A;
  Vec b,x;
  PC pc;
  KSP ksp;
  PetscInt mn,mi,mj,mnd,mid,mjd,ista,iend;
  PetscErrorCode ierr;
  std::fstream fid;
  void *ctx=0;
  float xmax,ymax,zmax,rdx,rdy,rdz;
  double dsend[3],drecv[3];

  PetscFunctionBegin;
  ierr = PetscInitialize(&argc, &argv, (char *) 0, (char *) 0);CHKERRQ(ierr);
  ncheck = 0;
  umem = 0;
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  if( myrank == 0 ) std::cout << "start" << std::endl;
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
  xq = new float [npmax];
  yq = new float [npmax];
  zq = new float [npmax];
  vq = new float [npmax];
  face = new int [3*npmax];
  idx = new int [npmax];
  vert = new float [3*npmax];
  are = new float [npmax];
  si = new double [npmax];
  un = new double [npmax];
  nfn = new int [npmax];
  na  = new int [npmax];
  nb  = new int [npmax];
  nc  = new int [npmax];
  nd  = new int [npmax];
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
  mem = npmax*26*4+npmax*3*8+100*8;
  memoryuse();

  // Read tripeptide data
  FILE* infile;
  char *dump;
  char fname[128],line[128];
  char fbase[128]="../../examples/biodata/lysozyme/copies_1000/lysozyme_1";
  char fxyzq[128]="../../examples/biodata/lysozyme/copies_1000/lysozyme";
  char fmove[128]="../../examples/biodata/lysozyme/copies_1000/gridMoves_1000_800.txt";
  int nmove;
  float bx,by,bz,dx,dy,dz,vx,vy,vz;
  double xt,yt,zt,th,ph;
  double xqd,yqd,zqd,vqd;
  sprintf(fname, "%s.vert", fbase);
  infile = fopen(fname, "r");
  nwn = 0;
  while (fgets(line, 128, infile)) nwn++;
  rewind(infile);
  for(PetscInt l = 0; l < nwn; ++l) {
    dump = fgets(line, 128, infile);
    sscanf(line, "%lf %lf %lf", &xt, &yt, &zt);
    vert[l*3+0] = xt; 
    vert[l*3+1] = yt;
    vert[l*3+2] = zt;
  }
  fclose(infile);
  sprintf(fname, "%s.face", fbase);
  infile = fopen(fname, "r");
  mi = 0;
  while (fgets(line, 128, infile)) mi++;
  rewind(infile);
  nwe = 0;
  for(PetscInt l = 0; l < mi; ++l) {
    dump = fgets(line, 128, infile);
    sscanf(line, "%u %u %u", &face[l*3+0], &face[l*3+1], &face[l*3+2]);
    bx = vert[3*face[3*l+2]-3]-vert[3*face[3*l+0]-3];
    dx = vert[3*face[3*l+1]-3]-vert[3*face[3*l+0]-3];
    by = vert[3*face[3*l+2]-2]-vert[3*face[3*l+0]-2];
    dy = vert[3*face[3*l+1]-2]-vert[3*face[3*l+0]-2];
    bz = vert[3*face[3*l+2]-1]-vert[3*face[3*l+0]-1];
    dz = vert[3*face[3*l+1]-1]-vert[3*face[3*l+0]-1];
    vx = by*dz-bz*dy;
    vy = bz*dx-bx*dz;
    vz = bx*dy-by*dx;
    if( 0.5*sqrt(vx*vx+vy*vy+vz*vz) > 1e-5 ) {
      face[3*nwe+0] = face[3*l+0];
      face[3*nwe+1] = face[3*l+1];
      face[3*nwe+2] = face[3*l+2];
      nwe++;
    }
  }
  fclose(infile);
  sprintf(fname, "%s.xyzq", fxyzq);
  infile = fopen(fname, "r");
  nch = 0;
  while (fgets(line, 128, infile)) nch++;
  rewind(infile);
  for(PetscInt l = 0; l < nch; ++l) {
    dump = fgets(line, 128, infile);
    sscanf(line, "%lf %lf %lf %lf", &xqd, &yqd, &zqd, &vqd);
    xq[l] = (float) xqd;
    yq[l] = (float) yqd;
    zq[l] = (float) zqd;
    vq[l] = (float) vqd;
  }
  fclose(infile);

  mn = nwn;
  mi = nwe;
  mj = nch;
  sprintf(fname, "%s", fmove);
  infile = fopen(fname, "r");
  nmove = 0;
  while (fgets(line, 128, infile)) nmove++;
  rewind(infile);
  nmove = 1;
  mpirange(0,nmove-1,ista,iend);
  for(PetscInt m = 0; m < nmove; ++m) {
    dump = fgets(line, 128, infile);
    sscanf(line, "%lf %lf %lf %lf %lf", &xt, &yt, &zt, &th, &ph);
    if( ista <= m && m <= iend ) {
      for(PetscInt l = 0; l < mi; ++l) {
        for(PetscInt k = 0; k < 3; ++k) face[3*nwe+k] = face[3*l+k]+nwn-mn;
        nwe++;
      }
      for(PetscInt l = 0; l < mn; ++l) {
        vert[3*nwn+0] = xt+cos(th)*cos(ph)*vert[3*l+0]-sin(th)*vert[3*l+1]+cos(th)*sin(ph)*vert[3*l+2];
        vert[3*nwn+1] = yt+sin(th)*cos(ph)*vert[3*l+0]+cos(th)*vert[3*l+1]+sin(th)*sin(ph)*vert[3*l+2];
        vert[3*nwn+2] = zt-sin(ph)*vert[3*l+0]+cos(ph)*vert[3*l+2];
        nwn++;
      }
      for(PetscInt l = 0; l < mj; ++l) {
        xq[nch] = xt+cos(th)*cos(ph)*xq[l]-sin(th)*yq[l]+cos(th)*sin(ph)*zq[l];
        yq[nch] = yt+sin(th)*cos(ph)*xq[l]+cos(th)*yq[l]+sin(th)*sin(ph)*zq[l];
        zq[nch] = zt-sin(ph)*xq[l]+cos(ph)*zq[l];
        vq[nch] = vq[l];
        nch++;
      }
    }
  }
  fclose(infile);

  mnd = 0;
  for(PetscInt l = mn; l < nwn; ++l) {
    for(PetscInt k = 0; k < 3; ++k) vert[3*mnd+k] = vert[3*l+k];
    mnd++;
  }
  mid = 0;
  for(PetscInt l = mi; l < nwe; ++l) {
    for(PetscInt k = 0; k < 3; ++k) face[3*mid+k] = face[3*l+k];
    mid++;
  }
  mjd = 0;
  for(PetscInt l = mj; l < nch; ++l) {
    xq[mjd] = xq[l];
    yq[mjd] = yq[l];
    zq[mjd] = zq[l];
    vq[mjd] = vq[l];
    mjd++;
  }
  nwn = mnd;
  nwe = mid;
  nch = mjd;

  mi = nwe;
  for(PetscInt l = 0; l < mi; ++l) {
    xi[l] = 0;
    yi[l] = 0;
    zi[l] = 0;
    for(PetscInt v = 0; v < 3; ++v) {
      xi[l] += vert[3*face[3*l+v]-3]/3;
      yi[l] += vert[3*face[3*l+v]-2]/3;
      zi[l] += vert[3*face[3*l+v]-1]/3;
    }
    bx = vert[3*face[3*l+2]-3]-vert[3*face[3*l+0]-3];
    dx = vert[3*face[3*l+1]-3]-vert[3*face[3*l+0]-3];
    by = vert[3*face[3*l+2]-2]-vert[3*face[3*l+0]-2];
    dy = vert[3*face[3*l+1]-2]-vert[3*face[3*l+0]-2];
    bz = vert[3*face[3*l+2]-1]-vert[3*face[3*l+0]-1];
    dz = vert[3*face[3*l+1]-1]-vert[3*face[3*l+0]-1];
    vx = by*dz-bz*dy;
    vy = bz*dx-bx*dz;
    vz = bx*dy-by*dx;
    are[l] = 0.5*sqrt(vx*vx+vy*vy+vz*vz);
    gxi[l] = vx/are[l]/2;
    gyi[l] = vy/are[l]/2;
    gzi[l] = vz/are[l]/2;
  }
  mj = nch;
  for(PetscInt l = 0; l < mj; ++l) {
    xj[l] = xq[l];
    yj[l] = yq[l];
    zj[l] = zq[l];
    vj[l] = vq[l];
  }
  for(PetscInt l = 0; l < 100; l++) tfmm[l] = 0;
  if(myrank==0) std::cout << "ni : " << mi << " nj : " << mj << std::endl;
#if 1
  if(myrank==0) std::cout << "direct" << std::endl;
  alloc(0);
  for(PetscInt l = 0; l < mi; ++l) vo[l] = 0;
  double tic,toc,tt=0;
  for(PetscInt m = 0; m < nprocs; ++m) {
    mpishift(mj,mjd,xj);
    mpishift(mj,mjd,yj);
    mpishift(mj,mjd,zj);
    mpishift(mj,mjd,vj);
    mj = mjd;
    tic = get_time();
    Gni_dir(0,mi-1,0,mj-1);
    toc = get_time();
    tt += toc-tic;
    for(PetscInt l = 0; l < mi; ++l) vo[l] += vi[l];
  }
  if(myrank==0) std::cout << "time : " << tt << std::endl;
  for(PetscInt l = 0; l < mi; ++l) vi[l] = vo[l];
  dealloc(0);
#else
  if(myrank==0) std::cout << "FMM" << std::endl;
  xmin = 1e6;
  xmax = -1e6;
  ymin = 1e6;
  ymax = -1e6;
  zmin = 1e6;
  zmax = -1e6; 
  for(PetscInt l = 0; l < mi; ++l) {
    xmin = std::min(xi[l],xmin);
    xmax = std::max(xi[l],xmax);
    ymin = std::min(yi[l],ymin);
    ymax = std::max(yi[l],ymax);
    zmin = std::min(zi[l],zmin);
    zmax = std::max(zi[l],zmax);
  } 
  for(PetscInt l = 0; l < mj; ++l) {
    xmin = std::min(xj[l],xmin);
    xmax = std::max(xj[l],xmax);
    ymin = std::min(yj[l],ymin);
    ymax = std::max(yj[l],ymax);
    zmin = std::min(zj[l],zmin);
    zmax = std::max(zj[l],zmax);
  } 
  dsend[0]=xmin;
  dsend[1]=ymin;
  dsend[2]=zmin;
  MPI_Allreduce(dsend,drecv,3,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  xmin = (float) drecv[0];
  ymin = (float) drecv[1];
  zmin = (float) drecv[2];
  dsend[0]=xmax;
  dsend[1]=ymax;
  dsend[2]=zmax;
  MPI_Allreduce(dsend,drecv,3,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  xmax = drecv[0];
  ymax = drecv[1]; 
  zmax = drecv[2];
  rdx = (xmax-xmin)*(1+1e-5);
  rdy = (ymax-ymin)*(1+1e-5);
  rdz = (zmax-zmin)*(1+1e-5);
  rd = std::max(rdx,rdy);
  rd = std::max(rd,rdz);
  nlevel(mj);
  boxpart(mi,mj);
  sorti(mi); 
  sortj(mj); 
  boxallocate(0,mi,0,mj);
  alloc(1);
  fmm(mi,mj,1);
  unsorti(mi); 
  unsortj(mj);
//  dealloc(1);
#endif
  MPI_Gather(tfmm,100,MPI_DOUBLE,tall,100,MPI_DOUBLE,0,MPI_COMM_WORLD);
  if( myrank == 0 ){
    fid.open("time.dat",std::ios::out | std::ios::binary);
    fid.write((char *)(&nprocs),sizeof(int));
    for(PetscInt j=0; j<nprocs; j++ ) {
      for(PetscInt i=0; i<100; i++ ) {
        fid.write((char *)(&tall[j][i]),sizeof(double));
      }
    }
    fid.close();
  }
  double sums = 0;
  for(PetscInt l=10; l < 23; ++l) sums += tfmm[l];
//  if(myrank==0) std::cout << "time : " << sums << std::endl;

  for(PetscInt l = 0; l < mi; ++l) {
    un[l] = -4*pi*are[l]*vi[l]/e1;
  }
  mj = mi;

// BEM-FMM Start here
  xmin = 1e6;
  xmax = -1e6;
  ymin = 1e6;
  ymax = -1e6;
  zmin = 1e6;
  zmax = -1e6;
  for(PetscInt l = 0; l < mi; ++l) {
    xmin = std::min(xi[l],xmin);
    xmax = std::max(xi[l],xmax);
    ymin = std::min(yi[l],ymin);
    ymax = std::max(yi[l],ymax);
    zmin = std::min(zi[l],zmin);
    zmax = std::max(zi[l],zmax);
  }
  dsend[0]=xmin;
  dsend[1]=ymin;
  dsend[2]=zmin;
  MPI_Allreduce(dsend,drecv,3,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  xmin = (float) drecv[0];
  ymin = (float) drecv[1];
  zmin = (float) drecv[2];
  dsend[0]=xmax;
  dsend[1]=ymax;
  dsend[2]=zmax;
  MPI_Allreduce(dsend,drecv,3,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  xmax = drecv[0];
  ymax = drecv[1];
  zmax = drecv[2];
  rdx = (xmax-xmin)*(1+1e-5);
  rdy = (ymax-ymin)*(1+1e-5);
  rdz = (zmax-zmin)*(1+1e-5);
  rd = std::max(rdx,rdy);
  rd = std::max(rd,rdz);
  for(PetscInt m = 0; m < mj; ++m) {
    xj[m] = xi[m];
    yj[m] = yi[m];
    zj[m] = zi[m];
    gxj[m] = (float) are[m];
    gyj[m] = 0;
    gzj[m] = (float) un[m];
  }
  nlevel(mj);
  boxpart(mi,mj);
  sorti(mi);
  sortj(mj);
  boxallocate(0,mi,0,mj);
  alloc(1);
  for(PetscInt m = 0; m < mj; ++m) {
    are[m] = gxj[m];
    si[m] = gyj[m];
    un[m] = gzj[m];
  }

  ierr = MatCreateShell(PETSC_COMM_WORLD,mi,mj,PETSC_DETERMINE,PETSC_DETERMINE,ctx,&A);CHKERRQ(ierr);
  ierr = MatShellSetOperation(A,MATOP_MULT, (void (*)(void)) mymatmult);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
  ierr = VecSetSizes(x,mi,PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&b);CHKERRQ(ierr);
  ierr = VecGetOwnershipRange(x,&ista,&iend);CHKERRQ(ierr);
  for(PetscInt l = 0; l < mi; ++l) idx[l] = l+ista;
  ierr = VecSetValues(x,mi,idx,si,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(x);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(x);CHKERRQ(ierr);
  ierr = VecSetValues(b,mi,idx,un,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCSHELL);CHKERRQ(ierr);
  ierr = PCShellSetApply(pc,mypcshell);CHKERRQ(ierr);
  ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
  ierr = VecGetValues(x,mi,idx,si);CHKERRQ(ierr);
  ierr = KSPDestroy(ksp);CHKERRQ(ierr);
  ierr = VecDestroy(x);CHKERRQ(ierr);
  ierr = VecDestroy(b);CHKERRQ(ierr);
  ierr = MatDestroy(A);CHKERRQ(ierr);

  for(PetscInt m = 0; m < mj; ++m) {
    gxj[m] = (float) are[m];
    gyj[m] = (float) si[m];
    gzj[m] = (float) un[m];
  }
  unsorti(mi);
  unsortj(mj);
  dealloc(1);
  for(PetscInt m = 0; m < mj; ++m) {
    are[m] = gxj[m];
    si[m] = gyj[m];
    un[m] = gzj[m];
  }

// End here

  for(PetscInt m = 0; m < mj; ++m) {
    xj[m] = xi[m];
    yj[m] = yi[m];
    zj[m] = zi[m];
    vj[m] = si[m]*are[m]/e1;
  }
  mi = nch;
  for(PetscInt l = 0; l < mi; ++l) {
    xi[l] = xq[l];
    yi[l] = yq[l];
    zi[l] = zq[l];
  }
#if 1
  alloc(0);
  for(PetscInt l = 0; l < mi; ++l) vo[l] = 0;
  for(PetscInt m = 0; m < nprocs; ++m) {
    mpishift(mj,mjd,xj);
    mpishift(mj,mjd,yj);
    mpishift(mj,mjd,zj);
    mpishift(mj,mjd,vj);
    mj = mjd;
    G_dir(0,mi-1,0,mj-1);
    for(PetscInt l = 0; l < mi; ++l) vo[l] += vi[l];
  }
  for(PetscInt l = 0; l < mi; ++l) vi[l] = vo[l];
  dealloc(0);
#else
  xmin = 1e6;
  xmax = -1e6;
  ymin = 1e6;
  ymax = -1e6;
  zmin = 1e6;
  zmax = -1e6;
  for(PetscInt l = 0; l < mi; ++l) {
    xmin = std::min(xi[l],xmin);
    xmax = std::max(xi[l],xmax);
    ymin = std::min(yi[l],ymin);
    ymax = std::max(yi[l],ymax);
    zmin = std::min(zi[l],zmin);
    zmax = std::max(zi[l],zmax);
  }
  for(PetscInt l = 0; l < mj; ++l) {
    xmin = std::min(xj[l],xmin);
    xmax = std::max(xj[l],xmax);
    ymin = std::min(yj[l],ymin);
    ymax = std::max(yj[l],ymax);
    zmin = std::min(zj[l],zmin);
    zmax = std::max(zj[l],zmax);
  } 
  dsend[0]=xmin;
  dsend[1]=ymin;
  dsend[2]=zmin;
  MPI_Allreduce(dsend,drecv,3,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  xmin = (float) drecv[0];
  ymin = (float) drecv[1];
  zmin = (float) drecv[2];
  dsend[0]=xmax;
  dsend[1]=ymax;
  dsend[2]=zmax;
  MPI_Allreduce(dsend,drecv,3,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  xmax = drecv[0];
  ymax = drecv[1];
  zmax = drecv[2];
  rdx = (xmax-xmin)*(1+1e-5);
  rdy = (ymax-ymin)*(1+1e-5);
  rdz = (zmax-zmin)*(1+1e-5);
  rd = std::max(rdx,rdy);
  rd = std::max(rd,rdz);
  nlevel(mj);
  boxpart(mi,mj);
  sorti(mi); 
  sortj(mj); 
  boxallocate(0,mi,0,mj);
  alloc(1);
  fmm(mi,mj,0);
  unsorti(mi);
  unsortj(mj);
  dealloc(1);
#endif
  for(PetscInt l = 0; l < mi; ++l) vi[l] *= 0.5*332.112*vq[l];

  sprintf(fname, "bio%d.dat", myrank);
  fid.open(fname,std::ios::out | std::ios::binary);

  fid.write((char *)(&nprocs),sizeof(int));
  fid.write((char *)(&nwe),sizeof(int));
  fid.write((char *)(&nwn),sizeof(int));
  fid.write((char *)(&nch),sizeof(int));
  for(PetscInt l = 0; l < 3*nwn; ++l) fid.write((char *)(&vert[l]),sizeof(double));
  for(PetscInt l = 0; l < 3*nwe; ++l) fid.write((char *)(&face[l]),sizeof(int));
  for(PetscInt l = 0; l < nwe; ++l) fid.write((char *)(&si[l]),sizeof(double));
  for(PetscInt l = 0; l < nch; ++l) fid.write((char *)(&xq[l]),sizeof(float));
  for(PetscInt l = 0; l < nch; ++l) fid.write((char *)(&yq[l]),sizeof(float));
  for(PetscInt l = 0; l < nch; ++l) fid.write((char *)(&zq[l]),sizeof(float));
  for(PetscInt l = 0; l < nch; ++l) fid.write((char *)(&vi[l]),sizeof(float));

  fid.close();
  ierr = PetscFinalize();CHKERRQ(ierr);

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
  delete[] xq;
  delete[] yq;
  delete[] zq;
  delete[] vq;
  delete[] face;
  delete[] idx;
  delete[] vert;
  delete[] are;
  delete[] si;
  delete[] un;
  delete[] nfn;
  delete[] na;
  delete[] nb;
  delete[] nc;
//  delete[] nd;
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

  mem = npmax*26*4+npmax*3*8+100*8;
  memoryfree();

  PetscFunctionReturn(0);
}
