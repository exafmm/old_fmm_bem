#define MAIN
#include "petscsnes.h"
#include "parameters.h"
#include "../include/constants.h"
#include "../include/arrays.h"
#undef MAIN

extern void memoryuse();
extern void memoryfree();
extern void geometry();
extern void gauss_trgl();
extern void elm_geom();
extern void nlevel(int);
extern void boxallocate(int, int, int, int);
extern void alloc(int);
extern void boxpart(int, int);
extern void sorti(int&);
extern void sortj(int&);
extern void unsorti(int&);
extern void unsortj(int&);
extern void dealloc(int);

#include "fmm_matmult.h"

int main(int argc, char *argv[]) {
  Mat A;
  Vec b,x;
  KSP ksp;
  PC pc;
  PetscInt i,ni,nj,nen,j=0,ij,ii;
  PetscScalar xd,yd,zd,qint,ph[6];
  PetscErrorCode ierr;
  void *ctx=0;
  std::fstream fid;

  PetscFunctionBegin;
  ierr = PetscInitialize(&argc, &argv, (char *) 0, (char *) 0);CHKERRQ(ierr);
  umem = 0;
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
  si = new double [npmax];
  un = new double [npmax];
  face = new int [npmax];
  vert = new float [npmax];
  xvert = new float [npmax][6];
  yvert = new float [npmax][6];
  zvert = new float [npmax][6];
  ptr = new float [npmax][3];
  ntr = new int [npmax][6];
  nne = new int [npmax][20];
  vna = new float [npmax][3];
  vta = new float [npmax][3];
  vsa = new float [npmax][3];
  p1 = new float [npmax][13];
  p2 = new float [npmax][13];
  p3 = new float [npmax][13];
  p4 = new float [npmax][13];
  p5 = new float [npmax][13];
  p6 = new float [npmax][13];
  hsg = new float [npmax][13];
  are = new float [npmax];
  crv = new float [npmax];
  xiq = new float [20];
  etq = new float [20];
  wq = new float [20];
  idx = new int [npmax];
  im = new int [npmax];
  nfn = new int [npmax];
  na  = new int [npmax];
  nb  = new int [npmax];
  nc  = new int [npmax];
  nd  = new int [npmax];
  nbi = new int [npmax];
  nbj = new int [npmax];
  sortd = new float [npmax];
  tfmm = new double [100];
  mem = npmax*180*4+100*8;
  memoryuse();

  geometry();
  gauss_trgl();
  elm_geom();

  std::cout << nwe << std::endl;
  for( i=0; i<nwn; i++ ) {
    xi[i] = ptr[i][0];
    yi[i] = ptr[i][1];
    zi[i] = ptr[i][2];
    gxi[i] = -vna[i][0];
    gyi[i] = -vna[i][1];
    gzi[i] = -vna[i][2];
  }
  ni = nwn;
  for( i=0; i<nch; i++ ) {
    xj[i] = xq[i];
    yj[i] = yq[i];
    zj[i] = zq[i];
    vj[i] = vq[i];
  }
  nj = nch;
  alloc(0);
  Gni_dir(0,ni-1,0,nj-1);
  dealloc(0);
  for( i=0; i<nwn; i++ ) {
    un[i] = -4*pi*vi[i]/e1;
    si[i] = 0;
  }

// Start here
  for( i=0; i<nwn; i++ ) idx[i] = i;
  ierr = MatCreateShell(PETSC_COMM_WORLD,nwn,nwn,PETSC_DECIDE,PETSC_DECIDE,ctx,&A);CHKERRQ(ierr);
  ierr = MatShellSetOperation(A,MATOP_MULT, (void (*)(void)) mymatmult);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
  ierr = VecSetSizes(x,PETSC_DECIDE,nwn);CHKERRQ(ierr);
  ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&b);CHKERRQ(ierr);
  ierr = VecSetValues(x,nwn,idx,si,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(x);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(x);CHKERRQ(ierr);
  ierr = VecSetValues(b,nwn,idx,un,INSERT_VALUES);CHKERRQ(ierr);
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
  ierr = VecGetValues(x,nwn,idx,si);CHKERRQ(ierr);
  ierr = KSPDestroy(ksp);CHKERRQ(ierr);
  ierr = VecDestroy(x);CHKERRQ(ierr);
  ierr = VecDestroy(b);CHKERRQ(ierr);
  ierr = MatDestroy(A);CHKERRQ(ierr);
// End here

  for( i=0; i<nch; i++ ) {
    xi[i] = xq[i];
    yi[i] = yq[i];
    zi[i] = zq[i];
  }
  ni = nch;
  nj = 0;
  // Run over source nodes
  for(PetscInt k = 0; k < nwe; ++k) {
    // Run over quadratures
    for(PetscInt j = 0; j < ngl; ++j) {
      xd    = 0;
      yd    = 0;
      zd    = 0;
      qint = 0;
      // This could be calculated on the fly from a reference element
      ph[0] = p1[k][j];
      ph[1] = p2[k][j];
      ph[2] = p3[k][j];
      ph[3] = p4[k][j];
      ph[4] = p5[k][j];
      ph[5] = p6[k][j];
      // Run over node contributions
      for(PetscInt i = 0; i < 6; ++i) {
        ii = ntr[k][i];
        // Quad point coordinates
        xd += ptr[ii][0]*ph[i];
        yd += ptr[ii][1]*ph[i];
        zd += ptr[ii][2]*ph[i];
        // Basis func value for active source point
        qint += ph[i]*si[ii]/e1;
      }
      // Collocation point
      xj[nj] = xd;
      yj[nj] = yd;
      zj[nj] = zd;
      vj[nj] = 0.5*hsg[k][j]*qint;
      nj++;
    }
  }

  alloc(0);
  G_dir(0,ni-1,0,nj-1);
  dealloc(0);
  for(PetscInt l = 0; l < ni; ++l) vi[l] *= 0.5*332.112*vq[l];
  float sums = 0;
  for(PetscInt l = 0; l < ni; ++l) sums += vi[l];
  std::cout << sums << std::endl;

  fid.open("bio.dat",std::ios::out | std::ios::binary);

  nen = 6;
  fid.write((char *)(&nwe),sizeof(int));
  fid.write((char *)(&nwn),sizeof(int));
  fid.write((char *)(&nen),sizeof(int));
  fid.write((char *)(&nch),sizeof(int));
  for( i=0; i<nwn; i++ ) fid.write((char *)(&ptr[i][0]),sizeof(float));
  for( i=0; i<nwn; i++ ) fid.write((char *)(&ptr[i][1]),sizeof(float));
  for( i=0; i<nwn; i++ ) fid.write((char *)(&ptr[i][2]),sizeof(float));
  for( i=0; i<nwn; i++ ) fid.write((char *)(&vna[i][0]),sizeof(double));
  for( i=0; i<nwn; i++ ) fid.write((char *)(&vna[i][1]),sizeof(double));
  for( i=0; i<nwn; i++ ) fid.write((char *)(&vna[i][2]),sizeof(double));
  for( i=0; i<nwn; i++ ) fid.write((char *)(&vta[i][0]),sizeof(double));
  for( i=0; i<nwn; i++ ) fid.write((char *)(&vta[i][1]),sizeof(double));
  for( i=0; i<nwn; i++ ) fid.write((char *)(&vta[i][2]),sizeof(double));
  for( i=0; i<nwn; i++ ) fid.write((char *)(&vsa[i][0]),sizeof(double));
  for( i=0; i<nwn; i++ ) fid.write((char *)(&vsa[i][1]),sizeof(double));
  for( i=0; i<nwn; i++ ) fid.write((char *)(&vsa[i][2]),sizeof(double));
  for( i=0; i<nch; i++ ) fid.write((char *)(&xq[i]),sizeof(float));
  for( i=0; i<nch; i++ ) fid.write((char *)(&yq[i]),sizeof(float));
  for( i=0; i<nch; i++ ) fid.write((char *)(&zq[i]),sizeof(float));
  for( i=0; i<nch; i++ ) fid.write((char *)(&vi[i]),sizeof(float));
  for( i=0; i<nwe; i++ ) {
    for( j=0; j<6; j++ ) {
      ij = i*6+j;
      idx[ij] = ntr[i][j];
    }
  }
  nwn = nwe*6;
  for( i=0; i<nwn; i++ ) fid.write((char *)(&ptr[idx[i]][0]),sizeof(double));
  for( i=0; i<nwn; i++ ) fid.write((char *)(&ptr[idx[i]][1]),sizeof(double));
  for( i=0; i<nwn; i++ ) fid.write((char *)(&ptr[idx[i]][2]),sizeof(double));
  for( i=0; i<nwn; i++ ) fid.write((char *)(&si[idx[i]]),sizeof(double));

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
  delete[] si;
  delete[] un;
  delete[] face;
  delete[] vert;
  delete[] xvert;
  delete[] yvert;
  delete[] zvert;
  delete[] ptr;
  delete[] ntr;
  delete[] nne;
  delete[] vna;
  delete[] vta;
  delete[] vsa;
  delete[] p1;
  delete[] p2;
  delete[] p3;
  delete[] p4;
  delete[] p5;
  delete[] p6;
  delete[] hsg;
  delete[] are;
  delete[] crv;
  delete[] xiq;
  delete[] etq;
  delete[] wq;
  delete[] idx;
  delete[] im;
  delete[] nfn;
  delete[] na;
  delete[] nb;
  delete[] nc;
  delete[] nd;
  delete[] nbi;
  delete[] nbj;
  delete[] sortd;
  delete[] tfmm;
  mem = npmax*180*4+100*8;
  memoryfree();

  PetscFunctionReturn(0);
}
