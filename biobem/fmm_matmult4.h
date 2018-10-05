extern int *im;
extern double *tfmm;

extern void G_dir(int, int, int, int);
extern void Gni_dir(int, int, int, int);
extern void Gnj_dir(int, int, int, int);
extern void fmm(int, int, int);

PetscErrorCode mymatmult(Mat A,Vec X,Vec Y)
{
  PetscInt i0,i1,i2,i3,i4,i5,i6,itest;
  PetscInt ii=0,ij,ni,nj;
  PetscScalar x0,y0,z0,ptl,x,y,z,vnx,vny,vnz,qint,area;
  PetscScalar dxij,dyij,dzij,r,den,cf,Gn;
  PetscScalar ph[6];
  PetscScalar *xx,*yy;
  PetscErrorCode ierr;
  float xmax,ymax,zmax,rdx,rdy,rdz;

  PetscFunctionBegin;
  ierr = VecSet(Y, 0.0);CHKERRQ(ierr);
  ierr = VecGetArray(X, &xx);CHKERRQ(ierr);
  ierr = VecGetArray(Y, &yy);CHKERRQ(ierr);

  for(PetscInt i = 0; i < 100; i++) tfmm[i] = 0;
  ncheck = 0;
  umem = 0;
  nprocs = 1;
  myrank = 0;

  xmin = 1e6;
  xmax = -1e6;
  ymin = 1e6;
  ymax = -1e6;
  zmin = 1e6;
  zmax = -1e6;
  // Prepare target data
  for(PetscInt l = 0; l < nwn; ++l) {
    xi[l] = ptr[l][0];
    yi[l] = ptr[l][1];
    zi[l] = ptr[l][2];
    xmin = std::min(xi[l],xmin);
    xmax = std::max(xi[l],xmax);
    ymin = std::min(yi[l],ymin);
    ymax = std::max(yi[l],ymax);
    zmin = std::min(zi[l],zmin);
    zmax = std::max(zi[l],zmax);
  }
  ni = nwn;

  for(PetscInt i = 0; i < nwn; ++i) im[i] = 0;
  // Run over source nodes
  nj = 0;
  for(PetscInt m = 0; m < nwn; ++m) {
    // Impulse
    im[m] = 1;
    // Run over source elements
    for(PetscInt kk = 1; kk <= nne[m][0]; ++kk) {
      PetscInt k = nne[m][kk];
      // Run over quadratures
      for(PetscInt j = 0; j < ngl; ++j) {
        x    = 0;
        y    = 0;
        z    = 0;
        vnx  = 0;
        vny  = 0;
        vnz  = 0;
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
          x += ptr[ii][0]*ph[i];
          y += ptr[ii][1]*ph[i];
          z += ptr[ii][2]*ph[i];
          // Normal at quad point
          vnx += vna[ii][0]*ph[i];
          vny += vna[ii][1]*ph[i];
          vnz += vna[ii][2]*ph[i];
          // Basis func value for active source point
          qint += im[ii]*ph[i];
        }
        xj[nj] = x;
        yj[nj] = y;
        zj[nj] = z;
        gxj[nj] = vnx;
        gyj[nj] = vny;
        gzj[nj] = vnz;
        vj[nj] = 0.5*hsg[k][j]*qint*xx[m]/e1;
        xmin = std::min(xj[nj],xmin);
        xmax = std::max(xj[nj],xmax);
        ymin = std::min(yj[nj],ymin);
        ymax = std::max(yj[nj],ymax);
        zmin = std::min(zj[nj],zmin);
        zmax = std::max(zj[nj],zmax);
        nj++;
      }
    }
    im[m] = 0; // Reset impulse
  }
  rdx = (xmax-xmin)*(1+1e-5);
  rdy = (ymax-ymin)*(1+1e-5);
  rdz = (zmax-zmin)*(1+1e-5);
  rd = std::max(rdx,rdy);
  rd = std::max(rd,rdz);

// Stage 1. Run for far field
  alloc(0);
  Gnj_dir(0,ni-1,0,nj-1);
  dealloc(0);
//  fmm(ni,nj,1);
  for(PetscInt l = 0; l < ni; ++l) {
    yy[l] = vi[l];
  }

  // Run over target nodes
  for(PetscInt l = 0; l < nwn; ++l) {
    x0 = ptr[l][0];
    y0 = ptr[l][1];
    z0 = ptr[l][2];
    ptl = 0;
    // Run over source elements
    for(PetscInt k = 0; k < nwe; ++k) {
      i1 = ntr[k][0];
      i2 = ntr[k][1];
      i3 = ntr[k][2];
      i4 = ntr[k][3];
      i5 = ntr[k][4];
      i6 = ntr[k][5];
      // Run over quadratures
      for(PetscInt j = 0; j < ngl; ++j) {
        x    = 0;
        y    = 0;
        z    = 0;
        vnx  = 0;
        vny  = 0;
        vnz  = 0;
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
          x += ptr[ii][0]*ph[i];
          y += ptr[ii][1]*ph[i];
          z += ptr[ii][2]*ph[i];
          // Normal at quad point
          vnx += vna[ii][0]*ph[i];
          vny += vna[ii][1]*ph[i];
          vnz += vna[ii][2]*ph[i];
        }
        // Compute the Green's function derivative and apply the triangle quadrature
        dxij = x0-x;
        dyij = y0-y;
        dzij = z0-z;
        r = sqrt(dxij*dxij+dyij*dyij+dzij*dzij);
        // Quad weight
        cf = 0.5*hsg[k][j];
        Gn = (vnx*dxij+vny*dyij+vnz*dzij)/(4*pi*r*r*r);
        ptl -= cf*Gn;
      }
    }
    yy[l] += ptl*xx[l]/e1;
  }

  // Run over target nodes
  for(PetscInt l = 0; l < nwn; ++l) {
    yy[l] += (e1+e2)/(e1-e2)*xx[l]/e1;
  }

  ierr = VecRestoreArray(X,&xx);CHKERRQ(ierr);
  ierr = VecRestoreArray(Y,&yy);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
