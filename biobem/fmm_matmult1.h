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
  for(PetscInt i = 0; i < nwn; ++i) im[i] = 0;
  // Run over source nodes
  nj = 0;
  for(PetscInt m = 0; m < nwn; ++m) {
    std::cout << m << " " << nwn << std::endl;
    // Impulse
    im[m] = 1;
    // Run over target nodes
    for(PetscInt l = 0; l < nwn; ++l) {
      x0 = ptr[l][0];
      y0 = ptr[l][1];
      z0 = ptr[l][2];
      i0 = im[l];
      ptl = 0;
      // Run over source elements
      for(PetscInt k = 0; k < nwe; ++k) {
        i1 = ntr[k][0];
        i2 = ntr[k][1];
        i3 = ntr[k][2];
        i4 = ntr[k][3];
        i5 = ntr[k][4];
        i6 = ntr[k][5];
        // Does this source element contain the current source node?
        itest = im[i1]+im[i2]+im[i3]+im[i4]+im[i5]+im[i6]+i0;
        if( itest != 0 ) {
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
            ph[0] = p1[j][k];
            ph[1] = p2[j][k];
            ph[2] = p3[j][k];
            ph[3] = p4[j][k];
            ph[4] = p5[j][k];
            ph[5] = p6[j][k];
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

            // Compute the Green's function derivative and apply the triangle quadrature
            dxij = x0-x;
            dyij = y0-y;
            dzij = z0-z;
            r = sqrt(dxij*dxij+dyij*dyij+dzij*dzij);
            den = 4*pi*r*r*r;
            // Quad weight
            cf = 0.5*hsg[k][j];

            Gn = vnx*dxij+vny*dyij+vnz*dzij;
            // This is the Green function contribution at the quad point, and then subtracting out the singular part
            ptl += (qint - i0)*Gn/den*cf;
          }
        }
      }
      yy[l] += (ptl+(e1+e2)/2/(e1-e2)*i0)*xx[m]/e1;
    }
    im[m] = 0.0; // Reset impulse
  }
  for(PetscInt i = 0; i < nwn; i++) {
    yy[i] += (e1+e2)/2/(e1-e2)*xx[i]/e1;
  }

  ierr = VecRestoreArray(X,&xx);CHKERRQ(ierr);
  ierr = VecRestoreArray(Y,&yy);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
