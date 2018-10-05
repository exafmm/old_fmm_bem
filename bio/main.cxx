#define MAIN
#include "petscsnes.h"
#include "parameters.h"
#include "../include/constants.h"
#include "../include/arrays.h"
#undef MAIN
#define DIR 0

extern void memoryuse();
extern void memoryfree();
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
  PetscInt ni,nj;
  PetscScalar *x,*y;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = VecSet(Y,0.0);CHKERRQ(ierr);
  ierr = VecGetArray(X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(Y,&y);CHKERRQ(ierr);

  ni = nwe;
  nj = nwe;
  for(PetscInt m = 0; m < nj; ++m) vj[m] = x[m]*are[m]/e1;

#if DIR
  Gni_dir(0,ni-1,0,nj-1);
#else
  fmm(ni,nj,1);
#endif

  for(PetscInt l = 0; l < ni; ++l) y[l] = vi[l]+(e1+e2)/(2*e1*(e1-e2))*x[l];

  ierr = VecRestoreArray(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(Y,&y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode mypcshell(PC pc,Vec X,Vec Y)
{
  PetscScalar *x,*y;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = VecSet(Y,0.0);CHKERRQ(ierr);
  ierr = VecGetArray(X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(Y,&y);CHKERRQ(ierr);

  for(PetscInt l = 0; l < nwe; ++l) y[l] = x[l]/(e1+e2)*(2*e1*(e1-e2));
//  for(PetscInt l = 0; l < nwe; ++l) y[l] = x[l];

  ierr = VecRestoreArray(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(Y,&y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

int main(int argc, char *argv[]) {
  Mat A;
  Vec b,x;
  KSP ksp;
  PC pc;
  PetscInt nn,ni,nj,nnd,nid,njd;
  PetscErrorCode ierr;
  std::fstream fid,fid2;
  void *ctx=0;
  float xmax,ymax,zmax,rdx,rdy,rdz;
  double tic,toc,t[4];

  PetscFunctionBegin;
  ierr = PetscInitialize(&argc, &argv, (char *) 0, (char *) 0);CHKERRQ(ierr);
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
  xq = new float [npmax];
  yq = new float [npmax];
  zq = new float [npmax];
  vq = new float [npmax];
  face = new int [npmax];
  idx = new int [npmax];
  vert = new float [npmax];
  are = new float [npmax];
  si = new double [npmax];
  un = new double [npmax];
  nfn = new int [npmax];
  na  = new int [npmax];
  nb  = new int [npmax];
  nc  = new int [npmax];
  nd  = new int [npmax];
  nbi = new int [npmax];
  nbj = new int [npmax];
  sortd = new float [npmax];
  tfmm = new double [100];
  mem = npmax*26*4+npmax*3*8+100*8;
  memoryuse();

  tic = get_time();
  // Read tripeptide data
  FILE* infile;
  char *dump;
  char fname[128],line[128];
//  char fbase[128]="../../examples/biodata/1oit_CDK_inhibitor/complex_1";
//  char fxyzq[128]="../../examples/biodata/1oit_CDK_inhibitor/complex";
//  char fmove[128]="../../examples/biodata/1oit_CDK_inhibitor/gridMoves_1000_800.txt";
  char fbase[128]="../../examples/biodata/lysozyme/copies_1000/lysozyme_1";
  char fxyzq[128]="../../examples/biodata/lysozyme/copies_1000/lysozyme";
//  char fmove[128]="../../examples/biodata/lysozyme/copies_1000/gridMoves_1000_800.txt";
//  char fbase[128]="../../examples/biodata/1oit_CDK_inhibitor/ligand_1";
//  char fxyzq[128]="../../examples/biodata/1oit_CDK_inhibitor/ligand";
//  char fbase[128]="../bio/data/sphere";
//  char fxyzq[128]="../bio/data/sphere";
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
  ni = 0;
  while (fgets(line, 128, infile)) ni++;
  rewind(infile);
  nwe = 0;
  for(PetscInt l = 0; l < ni; ++l) {
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

#if 0
  nn = nwn;
  ni = nwe;
  nj = nch;
  sprintf(fname, "%s", fmove);
  infile = fopen(fname, "r");
  nmove = 0;
  while (fgets(line, 128, infile)) nmove++;
  rewind(infile);
  nmove = 1;
  for(PetscInt m = 0; m < nmove; ++m) {
    dump = fgets(line, 128, infile);
    sscanf(line, "%lf %lf %lf %lf %lf", &xt, &yt, &zt, &th, &ph);
    for(PetscInt l = 0; l < ni; ++l) {
      for(PetscInt k = 0; k < 3; ++k) face[3*nwe+k] = face[3*l+k]+nwn-nn;
      nwe++;
    }
    for(PetscInt l = 0; l < nn; ++l) {
      vert[3*nwn+0] = xt+cos(th)*cos(ph)*vert[3*l+0]-sin(th)*vert[3*l+1]+cos(th)*sin(ph)*vert[3*l+2];
      vert[3*nwn+1] = yt+sin(th)*cos(ph)*vert[3*l+0]+cos(th)*vert[3*l+1]+sin(th)*sin(ph)*vert[3*l+2];
      vert[3*nwn+2] = zt-sin(ph)*vert[3*l+0]+cos(ph)*vert[3*l+2];
      nwn++;
    }
    for(PetscInt l = 0; l < nj; ++l) {
      xq[nch] = xt+cos(th)*cos(ph)*xq[l]-sin(th)*yq[l]+cos(th)*sin(ph)*zq[l];
      yq[nch] = yt+sin(th)*cos(ph)*xq[l]+cos(th)*yq[l]+sin(th)*sin(ph)*zq[l];
      zq[nch] = zt-sin(ph)*xq[l]+cos(ph)*zq[l];
      vq[nch] = vq[l];
      nch++;
    }
  }
  fclose(infile);
  
  nnd = 0;
  for(PetscInt l = nn; l < nwn; ++l) {
    for(PetscInt k = 0; k < 3; ++k) vert[3*nnd+k] = vert[3*l+k];
    nnd++;
  }
  nid = 0;
  for(PetscInt l = ni; l < nwe; ++l) {
    for(PetscInt k = 0; k < 3; ++k) face[3*nid+k] = face[3*l+k];
    nid++;
  }
  njd = 0;
  for(PetscInt l = nj; l < nch; ++l) {
    xq[njd] = xq[l];
    yq[njd] = yq[l];
    zq[njd] = zq[l];
    vq[njd] = vq[l];
    njd++;
  }
  nwn = nnd;
  nwe = nid;
  nch = njd;
#endif

  ni = nwe;
  for(PetscInt l = 0; l < ni; ++l) {
    bx = vert[3*face[3*l+2]-3]-vert[3*face[3*l+0]-3];
    dx = vert[3*face[3*l+1]-3]-vert[3*face[3*l+0]-3];
    by = vert[3*face[3*l+2]-2]-vert[3*face[3*l+0]-2];
    dy = vert[3*face[3*l+1]-2]-vert[3*face[3*l+0]-2];
    bz = vert[3*face[3*l+2]-1]-vert[3*face[3*l+0]-1];
    dz = vert[3*face[3*l+1]-1]-vert[3*face[3*l+0]-1];
    vx = by*dz-bz*dy;
    vy = bz*dx-bx*dz;
    vz = bx*dy-by*dx;
    xi[l] = 0;
    yi[l] = 0;
    zi[l] = 0;
    for(PetscInt v = 0; v < 3; ++v) {
      xi[l] += vert[3*face[3*l+v]-3]/3;
      yi[l] += vert[3*face[3*l+v]-2]/3;
      zi[l] += vert[3*face[3*l+v]-1]/3;
    }
    vx = by*dz-bz*dy;
    vy = bz*dx-bx*dz;
    vz = bx*dy-by*dx;
    are[l] = 0.5*sqrt(vx*vx+vy*vy+vz*vz);
    gxi[l] = vx/are[l]/2;
    gyi[l] = vy/are[l]/2;
    gzi[l] = vz/are[l]/2;
  }
  nj = nch;
  for(PetscInt l = 0; l < nj; ++l) {
    xj[l] = xq[l];
    yj[l] = yq[l];
    zj[l] = zq[l];
    vj[l] = vq[l];
  }
  toc = get_time();
  t[0] = toc-tic;
  tic = get_time();
#if DIR
  std::cout << "direct" << std::endl;
  alloc(0);
  Gni_dir(0,ni-1,0,nj-1);
  dealloc(0);
#else
  std::cout << "FMM" << std::endl;
  xmin = 1e6;
  xmax = -1e6;
  ymin = 1e6;
  ymax = -1e6;
  zmin = 1e6;
  zmax = -1e6;
  for(PetscInt l = 0; l < ni; ++l) {
    xmin = std::min(xi[l],xmin);
    xmax = std::max(xi[l],xmax);
    ymin = std::min(yi[l],ymin);
    ymax = std::max(yi[l],ymax);
    zmin = std::min(zi[l],zmin);
    zmax = std::max(zi[l],zmax);
  } 
  for(PetscInt l = 0; l < nj; ++l) {
    xmin = std::min(xj[l],xmin);
    xmax = std::max(xj[l],xmax);
    ymin = std::min(yj[l],ymin);
    ymax = std::max(yj[l],ymax);
    zmin = std::min(zj[l],zmin);
    zmax = std::max(zj[l],zmax);
  }
  rdx = (xmax-xmin)*(1+1e-5);
  rdy = (ymax-ymin)*(1+1e-5);
  rdz = (zmax-zmin)*(1+1e-5);
  rd = std::max(rdx,rdy);
  rd = std::max(rd,rdz);
  nlevel(nj);
  boxpart(ni,nj);
  sorti(ni);
  sortj(nj);
  boxallocate(0,ni,0,nj);
  alloc(1);
  fmm(ni,nj,1);
  unsorti(ni);
  unsortj(nj);
  dealloc(1);
#endif

  toc = get_time();
  t[1] = toc-tic;
  tic = get_time();
  
  for(PetscInt l = 0; l < ni; ++l) {
    un[l] = -4*pi*vi[l]/e1;
  }
  nj = ni;

// BEM-FMM Start here
  for(PetscInt i = 0; i < 100; i++) tfmm[i] = 0;
  xmin = 1e6;
  xmax = -1e6;
  ymin = 1e6;
  ymax = -1e6;
  zmin = 1e6;
  zmax = -1e6;
  for(PetscInt l = 0; l < ni; ++l) {
    xmin = std::min(xi[l],xmin);
    xmax = std::max(xi[l],xmax);
    ymin = std::min(yi[l],ymin);
    ymax = std::max(yi[l],ymax);
    zmin = std::min(zi[l],zmin);
    zmax = std::max(zi[l],zmax);
  }
  rdx = (xmax-xmin)*(1+1e-5);
  rdy = (ymax-ymin)*(1+1e-5);
  rdz = (zmax-zmin)*(1+1e-5);
  rd = std::max(rdx,rdy);
  rd = std::max(rd,rdz);
  for(PetscInt m = 0; m < nj; ++m) {
    xj[m] = xi[m];
    yj[m] = yi[m];
    zj[m] = zi[m];
    gxj[m] = (float) are[m];
    gyj[m] = 0;
    gzj[m] = (float) un[m];
  }
  nlevel(nj);
  boxpart(ni,nj);
  sorti(ni);
  sortj(nj);
  boxallocate(0,ni,0,nj);
  alloc(1);
  for(PetscInt m = 0; m < nj; ++m) {
    are[m] = gxj[m];
    si[m] = gyj[m];
    un[m] = gzj[m];
  }
  printf("N = %d\n",ni);

  for(PetscInt l = 0; l < ni; ++l) idx[l] = l;
  ierr = MatCreateShell(PETSC_COMM_WORLD,ni,nj,PETSC_DETERMINE,PETSC_DETERMINE,ctx,&A);CHKERRQ(ierr);
  ierr = MatShellSetOperation(A,MATOP_MULT, (void (*)(void)) mymatmult);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
  ierr = VecSetSizes(x,ni,PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&b);CHKERRQ(ierr);
  ierr = VecSetValues(x,ni,idx,si,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(x);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(x);CHKERRQ(ierr);
  ierr = VecSetValues(b,ni,idx,un,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCSHELL);CHKERRQ(ierr);
  ierr = PCShellSetApply(pc,mypcshell);CHKERRQ(ierr);
  ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  toc = get_time();
  t[0] += toc-tic;
  tic = get_time();
  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
  toc = get_time();
  t[2] = toc-tic;
  tic = get_time();
  ierr = VecGetValues(x,ni,idx,si);CHKERRQ(ierr);
  ierr = KSPDestroy(ksp);CHKERRQ(ierr);
  ierr = VecDestroy(x);CHKERRQ(ierr);
  ierr = VecDestroy(b);CHKERRQ(ierr);
  ierr = MatDestroy(A);CHKERRQ(ierr);

  for(PetscInt m = 0; m < nj; ++m) {
    gxj[m] = (float) are[m];
    gyj[m] = (float) si[m];
    gzj[m] = (float) un[m];
  }
  unsorti(ni);
  unsortj(nj);
  dealloc(1);
  for(PetscInt m = 0; m < nj; ++m) {
    are[m] = gxj[m];
    si[m] = gyj[m];
    un[m] = gzj[m];
  }
  fid.open("time.dat",std::ios::out | std::ios::binary);
  fid.write((char *)(&ni),sizeof(int));
  for(PetscInt i=0; i<100; i++ ) {
    fid.write((char *)(&tfmm[i]),sizeof(double));
  }
  fid.close();
// End here

  for(PetscInt m = 0; m < nj; ++m) {
    xj[m] = xi[m];
    yj[m] = yi[m];
    zj[m] = zi[m];
    vj[m] = si[m]*are[m]/e1;
  }
  ni = nch;
  for(PetscInt l = 0; l < ni; ++l) {
    xi[l] = xq[l];
    yi[l] = yq[l];
    zi[l] = zq[l];
  }
  toc = get_time();
  t[0] += toc-tic;
  tic = get_time();
#if DIR
  alloc(0);
  G_dir(0,ni-1,0,nj-1);
  dealloc(0);
#else
  xmin = 1e6;
  xmax = -1e6;
  ymin = 1e6;
  ymax = -1e6;
  zmin = 1e6;
  zmax = -1e6;
  for(PetscInt l = 0; l < ni; ++l) {
    xmin = std::min(xi[l],xmin);
    xmax = std::max(xi[l],xmax);
    ymin = std::min(yi[l],ymin);
    ymax = std::max(yi[l],ymax);
    zmin = std::min(zi[l],zmin);
    zmax = std::max(zi[l],zmax);
  }
  for(PetscInt l = 0; l < nj; ++l) {
    xmin = std::min(xj[l],xmin);
    xmax = std::max(xj[l],xmax);
    ymin = std::min(yj[l],ymin);
    ymax = std::max(yj[l],ymax);
    zmin = std::min(zj[l],zmin);
    zmax = std::max(zj[l],zmax);
  }
  rdx = (xmax-xmin)*(1+1e-5);
  rdy = (ymax-ymin)*(1+1e-5);
  rdz = (zmax-zmin)*(1+1e-5);
  rd = std::max(rdx,rdy);
  rd = std::max(rd,rdz);
  nlevel(nj);
  boxpart(ni,nj);
  sorti(ni);
  sortj(nj);
  boxallocate(0,ni,0,nj);
  alloc(1);
  fmm(ni,nj,0);
  unsorti(ni);
  unsortj(nj);
  dealloc(1);
#endif

  toc = get_time();
  t[3] = toc-tic;
  tic = get_time();
  for(PetscInt l = 0; l < ni; ++l) vi[l] *= 0.5*332.112*vq[l];
  float sums = 0;
  for(PetscInt l = 0; l < ni; ++l) sums += vi[l];
  std::cout << sums << std::endl;

  fid.open("bio0.dat",std::ios::out | std::ios::binary);

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

  for(PetscInt l = 0; l < nwe; ++l) fid.write((char *)(&xj[l]),sizeof(float));
  for(PetscInt l = 0; l < nwe; ++l) fid.write((char *)(&yj[l]),sizeof(float));
  for(PetscInt l = 0; l < nwe; ++l) fid.write((char *)(&zj[l]),sizeof(float));
  for(PetscInt l = 0; l < nwe; ++l) fid.write((char *)(&gxi[l]),sizeof(float));
  for(PetscInt l = 0; l < nwe; ++l) fid.write((char *)(&gyi[l]),sizeof(float));
  for(PetscInt l = 0; l < nwe; ++l) fid.write((char *)(&gzi[l]),sizeof(float));

  toc = get_time();
  t[0] += toc-tic;
  tic = get_time();
  fid2.open("time.dat",std::ios::out | std::ios::binary);
  fid2.write((char *)(&nj),sizeof(int));
  for(PetscInt l = 0; l < 100; ++l) fid2.write((char *)(&tfmm[l]),sizeof(double));
  fid2.close();

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
  delete[] nd;
  delete[] nbi;
  delete[] nbj;
  delete[] sortd;
  delete[] tfmm;
  mem = npmax*26*4+npmax*3*8+100*8;
  memoryfree();
  PetscFunctionReturn(0);
}
