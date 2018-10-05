#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <cutil.h>
#include <multithreading.h>

#include "var.h"
#include "G_p2p.cu"
#include "Gni_p2p.cu"
#include "Gnj_p2p.cu"
#include "G_p2m.cu"
#include "Gn_p2m.cu"
#include "m2m.cu"
#include "m2l.cu"
#include "l2l.cu"
#include "G_l2p.cu"
#include "Gn_l2p.cu"
#include "kernel.cu"
#include "wrapper.cu"

void G_p2p_gpu();
void Gni_p2p_gpu();
void Gnj_p2p_gpu();
void G_p2m_gpu();
void Gn_p2m_gpu();
void m2m_gpu();
void m2l_gpu();
void l2l_gpu();
void G_l2p_gpu();
void Gn_l2p_gpu();

double get_gpu_time(void)
{
  struct timeval tv;
  struct timezone tz;
  if (is_set==1) cudaThreadSynchronize();
  gettimeofday(&tv, &tz);
  return ((double)(tv.tv_sec+tv.tv_usec*1.0e-6));
}

SS add_nitadori(SS a, float ys)
{
  float ws;
  SS b;
  b.hs=a.hs+ys;
  ws=b.hs-a.hs;
  b.ls=ys-ws;
  b.ls=a.ls+b.ls;
  return b;
}

extern void p2pgpu_(int nvecd[],double *op,
	float xi[],float yi[],float zi[],
	float gxi[],float gyi[],float gzi[],float vi[],
	float xj[],float yj[],float zj[],
	float gxj[],float gyj[],float gzj[],float vj[])
{
  int i,nn,ni,nk,nflop;
  idev=nvecd[0];
  iblok=nvecd[1];
  int mblok=nvecd[2],nj=nvecd[3],neq=nvecd[4],myrank=nvecd[5];
  double tic,toc,flops,t[10];

  for(i=0;i<10;i++) t[i]=0;
  tic=get_gpu_time();
  nn=iblok*mblok+10;
  ni=iblok*NBLOK0;
  nk=iblok*NBLOK0;
  ms=sizeof(float)*NSCAL;
  mn=sizeof(int)*ROUNDUP0(nn);
  mi=sizeof(float)*6*ROUNDUP0(ni);
  mj=sizeof(float)*7*ROUNDUP0(nj);
  mk=sizeof(float)*4*ROUNDUP0(nk);

  scald=(float *)malloc(ms);
  nvec=(int *)malloc(mn);
  ivec=(float *)malloc(mi);
  jvec=(float *)malloc(mj);
  kvec=(float *)malloc(mk);

  if (is_set==0) {
    CUDA_SAFE_CALL(cudaSetDevice(idev));
    CUDA_SAFE_CALL(cudaGetDevice(&jdev));
    assert(idev==jdev);
    for(i=0;i<10;i++) tgpu[i]=0;
    is_set=1;
  }
  toc=tic;
  tic=get_gpu_time();
  t[0]+=tic-toc;
  if (mn>mn_a) {
    if(mn_a!=0) CUDA_SAFE_CALL(cudaFree(nveg[idev]));
    CUDA_SAFE_CALL(cudaMalloc((void**)(nveg+idev),mn));
    mn_a=mn;
  }
  if (mi>mi_a) {
    if(mi_a!=0) CUDA_SAFE_CALL(cudaFree(iveg[idev]));
    CUDA_SAFE_CALL(cudaMalloc((void**)(iveg+idev),mi));
    mi_a=mi;
  }
  if (mj>mj_a) {
    if(mj_a!=0) CUDA_SAFE_CALL(cudaFree(jveg[idev]));
    CUDA_SAFE_CALL(cudaMalloc((void**)(jveg+idev),mj));
    mj_a=mj;
  }
  if (mk>mk_a) {
    if(mk_a!=0) CUDA_SAFE_CALL(cudaFree(kveg[idev]));
    CUDA_SAFE_CALL(cudaMalloc((void**)(kveg+idev),mk));
    mk_a=mk;
  }
  toc=tic;
  tic=get_gpu_time();
  t[1]+=tic-toc;

  scald[0]=0;
  scald[1]=0;
  scald[2]=0;
  scald[3]=0;
  scald[4]=0;
  scald[5]=0;
  for(i=0;i<nn;i++){
    nvec[i]=nvecd[i];
  }
  for(i=0;i<ni;i++){
    ivec[6*i+0]=xi[i];
    ivec[6*i+1]=yi[i];
    ivec[6*i+2]=zi[i];
    ivec[6*i+3]=gxi[i];
    ivec[6*i+4]=gyi[i];
    ivec[6*i+5]=gzi[i];
  }
  for(i=0;i<nj;i++){
    jvec[7*i+0]=xj[i];
    jvec[7*i+1]=yj[i];
    jvec[7*i+2]=zj[i];
    jvec[7*i+3]=gxj[i];
    jvec[7*i+4]=gyj[i];
    jvec[7*i+5]=gzj[i];
    jvec[7*i+6]=vj[i];
  }

  toc=tic;
  tic=get_gpu_time();
  t[0]+=tic-toc;
  CUDA_SAFE_CALL(cudaMemcpyToSymbol(scal,scald,ms));
  CUDA_SAFE_CALL(cudaMemcpy(nveg[idev],nvec,mn,cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(iveg[idev],ivec,mi,cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(jveg[idev],jvec,mj,cudaMemcpyHostToDevice));
  toc=tic;
  tic=get_gpu_time();
  t[2]+=tic-toc;

  switch(neq){
    case 0 : G_p2p_gpu(); nflop=64 ;break;
    case 1 : Gni_p2p_gpu(); nflop=64 ;break;
    case 2 : Gnj_p2p_gpu(); nflop=64 ;break;
  }

  toc=tic;
  tic=get_gpu_time();
  t[3]+=tic-toc;
  CUDA_SAFE_CALL(cudaMemcpy(kvec,kveg[idev],mk,cudaMemcpyDeviceToHost));
  CUT_THREADEND;
  toc=tic;
  tic=get_gpu_time();
  t[2]+=tic-toc;

  for(i=0;i<nk;i++){
     vi[i]=kvec[4*i+0];
     gxi[i]=kvec[4*i+1];
     gyi[i]=kvec[4*i+2];
     gzi[i]=kvec[4*i+3];
  }

  free(scald);
  free(nvec);
  free(ivec);
  free(jvec);
  free(kvec);

  toc=tic;
  tic=get_gpu_time();
  t[0]+=tic-toc;
  for(i=0;i<9;i++) t[9]+=t[i];
  for(i=0;i<10;i++) tgpu[i]+=t[i];
  flops=*op*((double)nflop)/t[9];
//  printf("[%2d:%d] p2p cudaMalloc : %f s\n",myrank,idev,t[1]);
//  printf("[%2d:%d] p2p cudaMemcpy : %f s\n",myrank,idev,t[2]);
//  printf("[%2d:%d] p2p cudaKernel : %f s\n",myrank,idev,t[3]);
//  printf("[%2d:%d] p2p other      : %f s\n",myrank,idev,t[0]);
//  printf("[%2d:%d] p2p flops      : %f G\n",myrank,idev,flops/1e9);
  tic=flops+myrank;
}

extern void p2mgpu_(int nvecd[],
        double *op,double *rb,
        double *xmin,double *ymin,double *zmin,
        float  xjg[],float  yjg[],float  zjg[],
        float gxjg[],float gyjg[],float gzjg[], float vjg[],
        float brex[],float bimx[],float fac[])
{
  int i,nn,ni,nk,nflop;
  double tic,toc,flops,t[10];

  for(i=0;i<10;i++) t[i]=0;
  tic=get_gpu_time();
  idev=nvecd[0];
  iblok=nvecd[1];
  int mblok=nvecd[2],nj=nvecd[3],neq=nvecd[4],myrank=nvecd[5],mp=nvecd[6];
  nn=iblok*mblok+10;
  ni=iblok*NBLOK1;
  nk=mp*mp;
  ms=sizeof(float)*NSCAL;
  mn=sizeof(int)*ROUNDUP1(nn);
  mi=sizeof(float)*2*ROUNDUP1(ni);
  mj=sizeof(float)*7*ROUNDUP1(nj);
  mk=sizeof(float)*ROUNDUP1(nk);

  scald=(float *)malloc(ms);
  nvec=(int *)malloc(mn);
  ivec=(float *)malloc(mi);
  jvec=(float *)malloc(mj);
  kvec=(float *)malloc(mk);

  if (is_set==0) {
    CUDA_SAFE_CALL(cudaSetDevice(idev));
    CUDA_SAFE_CALL(cudaGetDevice(&jdev));
    assert(idev==jdev);
    for(i=0;i<10;i++) tgpu[i]=0;
    is_set=1;
  }
  toc=tic;
  tic=get_gpu_time();
  t[0]+=tic-toc;
  if (mn>mn_a) {
    if(mn_a!=0) CUDA_SAFE_CALL(cudaFree(nveg[idev]));
    CUDA_SAFE_CALL(cudaMalloc((void**)(nveg+idev),mn));
    mn_a=mn;
  }
  if (mi>mi_a) {
    if(mi_a!=0) CUDA_SAFE_CALL(cudaFree(iveg[idev]));
    CUDA_SAFE_CALL(cudaMalloc((void**)(iveg+idev),mi));
    mi_a=mi;
  }
  if (mj>mj_a) {
    if(mj_a!=0) CUDA_SAFE_CALL(cudaFree(jveg[idev]));
    CUDA_SAFE_CALL(cudaMalloc((void**)(jveg+idev),mj));
    mj_a=mj;
  }
  if (mk>mk_a) {
    if(mk_a!=0) CUDA_SAFE_CALL(cudaFree(kveg[idev]));
    CUDA_SAFE_CALL(cudaMalloc((void**)(kveg+idev),mk));
    mk_a=mk;
  }
  toc=tic;
  tic=get_gpu_time();
  t[1]+=tic-toc;

  scald[0]=(float)*rb;
  scald[1]=(float)*xmin;
  scald[2]=(float)*ymin;
  scald[3]=(float)*zmin;
  scald[4]=0;
  scald[5]=0;
  for(i=0;i<nn;i++){
    nvec[i]=nvecd[i];
  }
  for(i=0;i<nj;i++){
    jvec[7*i+0]= xjg[i];
    jvec[7*i+1]= yjg[i];
    jvec[7*i+2]= zjg[i];
    jvec[7*i+3]=gxjg[i];
    jvec[7*i+4]=gyjg[i];
    jvec[7*i+5]=gzjg[i];
    jvec[7*i+6]= vjg[i];
  }
  for(i=0;i<nk;i++){
    kvec[i]=fac[i];
  }

  toc=tic;
  tic=get_gpu_time();
  t[0]+=tic-toc;
  CUDA_SAFE_CALL(cudaMemcpyToSymbol(scal,scald,ms));
  CUDA_SAFE_CALL(cudaMemcpy(nveg[idev],nvec,mn,cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(jveg[idev],jvec,mj,cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(kveg[idev],kvec,mk,cudaMemcpyHostToDevice));
  toc=tic;
  tic=get_gpu_time();
  t[2]+=tic-toc;

  switch(neq){
    case 0 : G_p2m_gpu(); nflop=291 ;break;
    case 1 : Gn_p2m_gpu(); nflop=568 ;break;
  }

  toc=tic;
  tic=get_gpu_time();
  t[3]+=tic-toc;
  CUDA_SAFE_CALL(cudaMemcpy(ivec,iveg[idev],mi,cudaMemcpyDeviceToHost));
  CUT_THREADEND;
  toc=tic;
  tic=get_gpu_time();
  t[2]+=tic-toc;

  for(i=0;i<ni;i++){
    brex[i]=ivec[2*i+0];
    bimx[i]=ivec[2*i+1];
  }

  free(scald);
  free(nvec);
  free(ivec);
  free(jvec);
  free(kvec);

  toc=tic;
  tic=get_gpu_time();
  t[0]+=tic-toc;
  for(i=0;i<9;i++) t[9]+=t[i];
  for(i=0;i<10;i++) tgpu[i]+=t[i];
  flops=*op*((double)nflop)/t[9];
//  printf("[%2d:%d] p2m cudaMalloc : %f s\n",myrank,idev,t[1]);
//  printf("[%2d:%d] p2m cudaMemcpy : %f s\n",myrank,idev,t[2]);
//  printf("[%2d:%d] p2m cudaKernel : %f s\n",myrank,idev,t[3]);
//  printf("[%2d:%d] p2m other      : %f s\n",myrank,idev,t[0]);
//  printf("[%2d:%d] p2m flops      : %f G\n",myrank,idev,flops/1e9);
  tic=flops+myrank;
}

extern void m2lgpu_(int nvecd[],double *op,double *rb,
	float arex[],float aimx[],float brex[],float bimx[],
	float ynmre[],float ynmim[],float dnmre[],float dnmim[])
{
  int i,nn,ni,nk,nl,nflop;
  double tic,toc,flops,t[10];

  for(i=0;i<10;i++) t[i]=0;
  tic=get_gpu_time();
  idev=nvecd[0];
  iblok=nvecd[1];
  int mblok=nvecd[2],nj=nvecd[3],neq=nvecd[4],myrank=nvecd[5],mp=nvecd[6],nrbm=nvecd[7];
  int mpdnm=(4*mp*mp*mp-mp)/3;
  nn=iblok*mblok+10;
  ni=iblok*NBLOK1;
  nk=mp*mp;
  nl=mpdnm*2*nrbm;
  ms=sizeof(float)*NSCAL;
  mn=sizeof(int)*ROUNDUP1(nn);
  mi=sizeof(float)*2*ROUNDUP1(ni);
  mj=sizeof(float)*2*ROUNDUP1(nj);
  mk=sizeof(float)*2*ROUNDUP1(nk);
  ml=sizeof(float)*2*ROUNDUP1(nl);

  scald=(float *)malloc(ms);
  nvec=(int *)malloc(mn);
  ivec=(float *)malloc(mi);
  jvec=(float *)malloc(mj);
  kvec=(float *)malloc(mk);
  lvec=(float *)malloc(ml);

  if (is_set==0) {
    CUDA_SAFE_CALL(cudaSetDevice(idev));
    CUDA_SAFE_CALL(cudaGetDevice(&jdev));
    assert(idev==jdev);
    for(i=0;i<10;i++) tgpu[i]=0;
    is_set=1;
  }
  toc=tic;
  tic=get_gpu_time();
  t[0]+=tic-toc;
  if (mn>mn_a) {
    if(mn_a!=0) CUDA_SAFE_CALL(cudaFree(nveg[idev]));
    CUDA_SAFE_CALL(cudaMalloc((void**)(nveg+idev),mn));
    mn_a=mn;
  }
  if (mi>mi_a) {
    if(mi_a!=0) CUDA_SAFE_CALL(cudaFree(iveg[idev]));
    CUDA_SAFE_CALL(cudaMalloc((void**)(iveg+idev),mi));
    mi_a=mi;
  }
  if (mj>mj_a) {
    if(mj_a!=0) CUDA_SAFE_CALL(cudaFree(jveg[idev]));
    CUDA_SAFE_CALL(cudaMalloc((void**)(jveg+idev),mj));
    mj_a=mj;
  }
  if (mk>mk_a) {
    if(mk_a!=0) CUDA_SAFE_CALL(cudaFree(kveg[idev]));
    CUDA_SAFE_CALL(cudaMalloc((void**)(kveg+idev),mk));
    mk_a=mk;
  }
  if (ml>ml_a) {
    if(ml_a!=0) CUDA_SAFE_CALL(cudaFree(lveg[idev]));
    CUDA_SAFE_CALL(cudaMalloc((void**)(lveg+idev),ml));
    ml_a=ml;
  }
  toc=tic;
  tic=get_gpu_time();
  t[1]+=tic-toc;

  scald[0]=(float)*rb;
  scald[1]=0;
  scald[2]=0;
  scald[3]=0;
  scald[4]=0;
  scald[5]=0;
  for(i=0;i<nn;i++){
    nvec[i]=nvecd[i];
  }
  for(i=0;i<nj;i++){
    jvec[2*i+0]=brex[i];
    jvec[2*i+1]=bimx[i];
  }
  for(i=0;i<nk;i++){
    kvec[2*i+0]=ynmre[i];
    kvec[2*i+1]=ynmim[i];
  }
  for(i=0;i<nl;i++){
    lvec[2*i+0]=dnmre[i];
    lvec[2*i+1]=dnmim[i];
  }

  toc=tic;
  tic=get_gpu_time();
  t[0]+=tic-toc;
  CUDA_SAFE_CALL(cudaMemcpyToSymbol(scal,scald,ms));
  CUDA_SAFE_CALL(cudaMemcpy(nveg[idev],nvec,mn,cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(jveg[idev],jvec,mj,cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(kveg[idev],kvec,mk,cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(lveg[idev],lvec,ml,cudaMemcpyHostToDevice));
  toc=tic;
  tic=get_gpu_time();
  t[2]+=tic-toc;

  switch(neq){
    case 0 : m2m_gpu(); nflop=48 ;break;
    case 1 : m2l_gpu(); nflop=48 ;break;
    case 2 : l2l_gpu(); nflop=48 ;break;
  }

  toc=tic;
  tic=get_gpu_time();
  t[3]+=tic-toc;
  CUDA_SAFE_CALL(cudaMemcpy(ivec,iveg[idev],mi,cudaMemcpyDeviceToHost));
  CUT_THREADEND;
  toc=tic;
  tic=get_gpu_time();
  t[2]+=tic-toc;

  for(i=0;i<ni;i++){
    arex[i]=ivec[2*i+0];
    aimx[i]=ivec[2*i+1];
  }

  free(scald);
  free(nvec);
  free(ivec);
  free(jvec);
  free(kvec);
  free(lvec);

  toc=tic;
  tic=get_gpu_time();
  t[0]+=tic-toc;
  for(i=0;i<9;i++) t[9]+=t[i];
  for(i=0;i<10;i++) tgpu[i]+=t[i];
  flops=*op*((double)nflop)/t[9];
//  printf("[%2d:%d] m2l cudaMalloc : %f s\n",myrank,idev,t[1]);
//  printf("[%2d:%d] m2l cudaMemcpy : %f s\n",myrank,idev,t[2]);
//  printf("[%2d:%d] m2l cudaKernel : %f s\n",myrank,idev,t[3]);
//  printf("[%2d:%d] m2l other      : %f s\n",myrank,idev,t[0]);
//  printf("[%2d:%d] m2l flops      : %f G\n",myrank,idev,flops/1e9);
  tic=flops+myrank;
}

extern void l2pgpu_(int nvecd[],
        double *op,double *rb,
        double *xmin,double *ymin,double *zmin,
        float  xig[],float  yig[],float  zig[],
        float gxig[],float gyig[],float gzig[],float vig[],
        float arex[],float aimx[],float fac[])
{
  int i,nn,ni,nk,nl,nflop;
  double tic,toc,flops,t[10];

  for(i=0;i<10;i++) t[i]=0;
  tic=get_gpu_time();
  idev=nvecd[0];
  iblok=nvecd[1];
  int mblok=nvecd[2],nj=nvecd[3],neq=nvecd[4],myrank=nvecd[5],mp=nvecd[6];
  nn=iblok*mblok+10;
  ni=iblok*NBLOK1;
  nk=iblok*NBLOK1;
  nl=mp*mp;
  ms=sizeof(float)*NSCAL;
  mn=sizeof(int)*ROUNDUP1(nn);
  mi=sizeof(float)*6*ROUNDUP1(ni);
  mj=sizeof(float)*2*ROUNDUP1(nj);
  mk=sizeof(float)*4*ROUNDUP1(nk);
  ml=sizeof(float)*ROUNDUP1(nl);

  scald=(float *)malloc(ms);
  nvec=(int *)malloc(mn);
  ivec=(float *)malloc(mi);
  jvec=(float *)malloc(mj);
  kvec=(float *)malloc(mk);
  lvec=(float *)malloc(ml);

  if (is_set==0) {
    CUDA_SAFE_CALL(cudaSetDevice(idev));
    CUDA_SAFE_CALL(cudaGetDevice(&jdev));
    assert(idev==jdev);
    for(i=0;i<10;i++) tgpu[i]=0;
    is_set=1;
  }
  toc=tic;
  tic=get_gpu_time();
  t[0]+=tic-toc;
  if (mn>mn_a) {
    if(mn_a!=0) CUDA_SAFE_CALL(cudaFree(nveg[idev]));
    CUDA_SAFE_CALL(cudaMalloc((void**)(nveg+idev),mn));
    mn_a=mn;   
  }
  if (mi>mi_a) {
    if(mi_a!=0) CUDA_SAFE_CALL(cudaFree(iveg[idev]));
    CUDA_SAFE_CALL(cudaMalloc((void**)(iveg+idev),mi));
    mi_a=mi;
  }
  if (mj>mj_a) {
    if(mj_a!=0) CUDA_SAFE_CALL(cudaFree(jveg[idev]));
    CUDA_SAFE_CALL(cudaMalloc((void**)(jveg+idev),mj));
    mj_a=mj;
  }
  if (mk>mk_a) {
    if(mk_a!=0) CUDA_SAFE_CALL(cudaFree(kveg[idev]));
    CUDA_SAFE_CALL(cudaMalloc((void**)(kveg+idev),mk));
    mk_a=mk;
  } 
  if (ml>ml_a) {
    if(ml_a!=0) CUDA_SAFE_CALL(cudaFree(lveg[idev]));
    CUDA_SAFE_CALL(cudaMalloc((void**)(lveg+idev),ml));
    ml_a=ml;
  }
  toc=tic;
  tic=get_gpu_time();
  t[1]+=tic-toc;

  scald[0]=(float)*rb;
  scald[1]=(float)*xmin;
  scald[2]=(float)*ymin;
  scald[3]=(float)*zmin;
  scald[4]=0;
  scald[5]=0;
  for(i=0;i<nn;i++){
    nvec[i]=nvecd[i];
  }
  for(i=0;i<ni;i++){
    ivec[6*i+0]= xig[i];
    ivec[6*i+1]= yig[i];
    ivec[6*i+2]= zig[i];
    ivec[6*i+3]=gxig[i];
    ivec[6*i+4]=gyig[i];
    ivec[6*i+5]=gzig[i];
  }
  for(i=0;i<nj;i++){
    jvec[2*i+0]=arex[i];
    jvec[2*i+1]=aimx[i];
  }
  for(i=0;i<nl;i++){
    lvec[i]=fac[i];
  }

  toc=tic;
  tic=get_gpu_time();
  t[0]+=tic-toc;
  CUDA_SAFE_CALL(cudaMemcpyToSymbol(scal,scald,ms));
  CUDA_SAFE_CALL(cudaMemcpy(nveg[idev],nvec,mn,cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(iveg[idev],ivec,mi,cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(jveg[idev],jvec,mj,cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(lveg[idev],lvec,ml,cudaMemcpyHostToDevice));
  toc=tic;
  tic=get_gpu_time();
  t[2]+=tic-toc;

  switch(neq){
    case 0 : G_l2p_gpu(); nflop=169 ;break;
    case 1 : Gn_l2p_gpu(); nflop=186 ;break;
  }

  toc=tic;
  tic=get_gpu_time();
  t[3]+=tic-toc;
  CUDA_SAFE_CALL(cudaMemcpy(kvec,kveg[idev],mk,cudaMemcpyDeviceToHost));
  CUT_THREADEND;
  toc=tic;
  tic=get_gpu_time();
  t[2]+=tic-toc;

  for(i=0;i<nk;i++){
     gxig[i]=kvec[4*i+0];
     gyig[i]=kvec[4*i+1];
     gzig[i]=kvec[4*i+2];
     vig[i]=kvec[4*i+3];
  }

  free(scald);
  free(nvec);
  free(ivec);
  free(jvec);
  free(kvec);
  free(lvec);

  toc=tic;
  tic=get_gpu_time();
  t[0]+=tic-toc;
  for(i=0;i<9;i++) t[9]+=t[i];
  for(i=0;i<10;i++) tgpu[i]+=t[i];
  flops=*op*((double)nflop)/t[9];
//  printf("[%2d:%d] l2p cudaMalloc : %f s\n",myrank,idev,t[1]);
//  printf("[%2d:%d] l2p cudaMemcpy : %f s\n",myrank,idev,t[2]);
//  printf("[%2d:%d] l2p cudaKernel : %f s\n",myrank,idev,t[3]);
//  printf("[%2d:%d] l2p other      : %f s\n",myrank,idev,t[0]);
//  printf("[%2d:%d] l2p flops      : %f G\n",myrank,idev,flops/1e9);
//  printf("[%2d:%d] l2p cudaMalloc : %f s\n",myrank,idev,tgpu[1]);
//  printf("[%2d:%d] l2p cudaMemcpy : %f s\n",myrank,idev,tgpu[2]);
//  printf("[%2d:%d] l2p cudaKernel : %f s\n",myrank,idev,tgpu[3]);
//  printf("[%2d:%d] l2p other      : %f s\n",myrank,idev,tgpu[0]);
  tic=flops+myrank;
}
