#include "../include/constants.h"

extern float *xq,*yq,*zq,*vq;
extern int *face;
extern float *vert;
extern float (*xvert)[6],(*yvert)[6],(*zvert)[6];
extern float (*ptr)[3];
extern int (*ntr)[6],(*nne)[20];

void geometry() {
  FILE* infile;
  char *dump;
  char fname[128],line[128];
//  char fbase[128]="../../examples/biodata/tripeptide/tripeptide_1";
//  char fxyzq[128]="../../examples/biodata/tripeptide/tripeptide";
//  char fbase[128]="../../examples/biodata/1oit_CDK_inhibitor/ligand_1";
//  char fxyzq[128]="../../examples/biodata/1oit_CDK_inhibitor/ligand";
//  char fbase[128]="../../examples/biodata/1oit_CDK_inhibitor/receptor_1";
//  char fxyzq[128]="../../examples/biodata/1oit_CDK_inhibitor/receptor";
  char fbase[128]="../bio/data/sphere";
  char fxyzq[128]="../bio/data/sphere";
  int i,nen,ndv,k,j,ic,ni;
  double xd,yd,zd,vd,bx,by,bz,dx,dy,dz,vx,vy,vz;

  // Read tripeptide data
  sprintf(fname, "%s.vert", fbase);
  infile = fopen(fname, "r");
  nwn = 0;
  while (fgets(line, 128, infile)) nwn++;
  rewind(infile);
  for(int i = 0; i < nwn; ++i) {
    dump = fgets(line, 128, infile);
    sscanf(line, "%lf %lf %lf", &xd, &yd, &zd);
    vert[i*3+0] = (float) xd;
    vert[i*3+1] = (float) yd;
    vert[i*3+2] = (float) zd;
  }
  fclose(infile);
  sprintf(fname, "%s.face", fbase);
  infile = fopen(fname, "r");
  ni = 0;
  while (fgets(line, 128, infile)) ni++;
  rewind(infile);
  nwe = 0;
  for(int i = 0; i < ni; ++i) {
    dump = fgets(line, 128, infile);
    sscanf(line, "%u %u %u", &face[i*3+0], &face[i*3+1], &face[i*3+2]);
    bx = vert[3*face[3*i+2]-3]-vert[3*face[3*i+0]-3];
    dx = vert[3*face[3*i+1]-3]-vert[3*face[3*i+0]-3];
    by = vert[3*face[3*i+2]-2]-vert[3*face[3*i+0]-2];
    dy = vert[3*face[3*i+1]-2]-vert[3*face[3*i+0]-2];
    bz = vert[3*face[3*i+2]-1]-vert[3*face[3*i+0]-1];
    dz = vert[3*face[3*i+1]-1]-vert[3*face[3*i+0]-1];
    vx = by*dz-bz*dy;
    vy = bz*dx-bx*dz;
    vz = bx*dy-by*dx;
    if( 0.5*sqrt(vx*vx+vy*vy+vz*vz) > 1e-5 ) {
      face[3*nwe+0] = face[3*i+0];
      face[3*nwe+1] = face[3*i+1];
      face[3*nwe+2] = face[3*i+2];
      nwe++;
    }
  }
  fclose(infile);
  sprintf(fname, "%s.xyzq", fxyzq);
  infile = fopen(fname, "r");
  nch = 0;
  while (fgets(line, 128, infile)) nch++;
  rewind(infile);
  for(int i = 0; i < nch; ++i) {
    dump = fgets(line, 128, infile);
    sscanf(line, "%lf %lf %lf %lf", &xd, &yd, &zd, &vd);
    xq[i] = (float) xd;
    yq[i] = (float) yd;
    zq[i] = (float) zd;
    vq[i] = (float) vd;
  }
  fclose(infile);
  for(int i = 0; i < nwe; ++i) {
    for(int v = 0; v < 3; ++v) {
      xvert[i][v] = vert[3*face[3*i+v]-3];
      yvert[i][v] = vert[3*face[3*i+v]-2];
      zvert[i][v] = vert[3*face[3*i+v]-1];
    }
  }

  nen = 6;
  ndv = 0;

  // Compute the midpoints
  for( i=0; i<nwe; i++ ) {
    xvert[i][3] = 0.5*(xvert[i][0]+xvert[i][1]);
    yvert[i][3] = 0.5*(yvert[i][0]+yvert[i][1]);
    zvert[i][3] = 0.5*(zvert[i][0]+zvert[i][1]);
    xvert[i][4] = 0.5*(xvert[i][1]+xvert[i][2]);
    yvert[i][4] = 0.5*(yvert[i][1]+yvert[i][2]);
    zvert[i][4] = 0.5*(zvert[i][1]+zvert[i][2]);
    xvert[i][5] = 0.5*(xvert[i][2]+xvert[i][0]);
    yvert[i][5] = 0.5*(yvert[i][2]+yvert[i][0]);
    zvert[i][5] = 0.5*(zvert[i][2]+zvert[i][0]);
  }

  // ptr[j][i]: Coordinate of the ith component for jth global node
  for(i = 0; i < nen; i++) {
    ptr[i][0] = xvert[0][i];
    ptr[i][1] = yvert[0][i];
    ptr[i][2] = zvert[0][i];
    ntr[0][i] = i;
  }
  nwn = nen;

  // ntr[j][i]: Global label on jth node on ith element
  for(i = 1; i < nwe; i++) {
    for(j = 0; j < nen; j++) {
      ic = 0;
      for(k = 0; k < nwn; k++) {
        if (std::abs(xvert[i][j]-ptr[k][0]) <= eps) {
          if (std::abs(yvert[i][j]-ptr[k][1]) <= eps) {
            if (std::abs(zvert[i][j]-ptr[k][2]) <= eps) {
              ic = 1;
              ntr[i][j] = k;
            }
          }
        }
      }
      if (ic == 0) {
        ptr[nwn][0] = xvert[i][j];
        ptr[nwn][1] = yvert[i][j];
        ptr[nwn][2] = zvert[i][j];
        ntr[i][j] = nwn;
        nwn++;
      }
    }
  }
  
  // nne[i][0] is the number of elements touching  global node i
  // nne[i][j] for j = 2,...,nne[i][0]+1 are the corresponding element labels
  for(i = 0; i < nwn; i++) {
    for(j = 0; j < 8; j++) {
      nne[i][j] = 0;
    }
  }
  for(i = 0; i < nwn; i++) {
    nne[i][0] = 0;
    ic = 0;
    for(j = 0; j < nwe; j++) {
      for(k = 0; k < nen; k++) {
        if (ntr[j][k] == i) {
          nne[i][0]++;
          ic++;
          nne[i][ic] = j;
        }
      }
    }
  }
}
