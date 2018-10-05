#define P2P_KERNEL(p2p_kernel,p2p_kernel_core) \
extern "C" __global__ \
void p2p_kernel(int *nvec,float *ivec,float *jvec,float *kvec)\
{\
  int bx=blockIdx.x;\
  int tx=threadIdx.x;\
  int mblok=nvec[2];\
  int ib,jb,jbase,jsize,jblok,nij;\
  int tx7,jj7;\
  int i,j,ij,jj;\
  float dxij,dyij,dzij,rij,rsij=0;\
  float eps=1e-6;\
  float pi14=0.25/M_PI;\
  float veci[6],veck[4];\
  __shared__ float vecj[NBLOK0*7];\
  rij=rsij;\
  ib=bx*NBLOK0+tx;\
  for(i=0;i<6;i++) veci[i]=ivec[6*ib+i];\
  for(i=0;i<4;i++) veck[i]=0.0f;\
  tx7=tx*7;\
  nij=nvec[bx*mblok+10];\
  for(ij=0;ij<nij;ij++){\
    jbase=nvec[bx*mblok+2*ij+11];\
    jsize=nvec[bx*mblok+2*ij+12];\
    jblok=(jsize+NBLOK0-1)/NBLOK0;\
    for(j=0;j<jblok-1;j++){\
      jb=j*NBLOK0+jbase+tx;\
      for(i=0;i<7;i++) vecj[tx7+i]=jvec[7*jb+i];\
      __syncthreads();\
      for(jj=jj7=0;jj<NBLOK0;jj++){\
        p2p_kernel_core;\
      }\
      __syncthreads();\
    }\
    jb=j*NBLOK0+jbase+tx;\
    for(i=0;i<7;i++) vecj[tx7+i]=jvec[7*jb+i];\
    jb=j*NBLOK1+jbase;\
    __syncthreads();\
    for(jj=jj7=0;jj<jsize-(j*NBLOK0);jj++){\
      p2p_kernel_core;\
    }\
    __syncthreads();\
  }\
  for(i=0;i<4;i++) kvec[4*ib+i]=veck[i];\
}

#define P2M_KERNEL(p2m_kernel,p2m_kernel_core) \
extern "C" __global__ \
void p2m_kernel(int *nvec,float *ivec,float *jvec,float *kvec)\
{\
  int bx=blockIdx.x;\
  int tx=threadIdx.x;\
  int i,j,k,m,n,ib,jb,jj,jj7,jbase,jsize,jblok,nm,nms;\
  int mblok=nvec[2],mp=nvec[6],nb,nc[3],nd;\
  float rb=scal[0],xmin=scal[1],ymin=scal[2],zmin=scal[3];\
  float xjc,yjc,zjc,xjjc,yjjc,zjjc,rh,al,be,eps=1e-6;\
  float xx,s2,p,pn,p1,p2,fact,ere,eim,rhm,rhn;\
  __shared__ int mg[NBLOK1],ng[NBLOK1];\
  __shared__ float veci[2*NBLOK1],vecj[7*NBLOK1],veck[MPMAX];\
  __shared__ float bnm[MPMAX];\
  ib=bx*NBLOK1+tx;\
  for(i=0;i<NBLOK1;i++){\
    ng[i]=0;\
    mg[i]=0;\
  }\
  for(n=0;n<mp;n++){\
    for(m=0;m<=n;m++){\
      nms=n*(n+1)/2+m;\
      ng[nms]=n;\
      mg[nms]=m;\
    }\
  }\
  jblok=(MPMAX+NBLOK1-1)/NBLOK1;\
  for(j=0;j<jblok-1;j++){\
    jb=j*NBLOK1+tx;\
    veck[jb]=kvec[jb];\
    veck[jb]=kvec[jb];\
    __syncthreads();\
  }\
  if(j*NBLOK1+tx<MPMAX){\
    jb=j*NBLOK1+tx;\
    veck[jb]=kvec[jb];\
    veck[jb]=kvec[jb];\
  }\
  __syncthreads();\
  for(i=0;i<2;i++) veci[2*tx+i]=0;\
  __syncthreads();\
  nb=nvec[bx*mblok+10];\
  jbase=nvec[bx*mblok+11];\
  jsize=nvec[bx*mblok+12];\
  for(i=0;i<3;i++) nc[i]=0;\
  k=0;\
  i=1;\
  while(nb!=0){\
    j=2-k;\
    nc[j]=nc[j]+nb%2*i;\
    nb=nb/2;\
    j=k+1;\
    k=j%3;\
    if(k==0) i=i*2;\
  }\
  nd=nc[0];\
  nc[0]=nc[1];\
  nc[1]=nc[2];\
  nc[2]=nd;\
  xjc=xmin+(nc[0]+0.5)*rb;\
  yjc=ymin+(nc[1]+0.5)*rb;\
  zjc=zmin+(nc[2]+0.5)*rb;\
  jblok=(jsize+NBLOK1-1)/NBLOK1;\
  for(j=0;j<jblok-1;j++){\
    jb=j*NBLOK1+jbase+tx;\
    for(i=0;i<7;i++) vecj[7*tx+i]=jvec[7*jb+i];\
    __syncthreads();\
    for(jj=jj7=0;jj<NBLOK1;jj++){\
      p2m_kernel_core;\
      __syncthreads();\
    }\
  }\
  jb=j*NBLOK1+jbase+tx;\
  for(i=0;i<7;i++) vecj[7*tx+i]=jvec[7*jb+i];\
  __syncthreads();\
  for(jj=jj7=0;jj<jsize-(j*NBLOK1);jj++){\
    p2m_kernel_core;\
    __syncthreads();\
  }\
  for(i=0;i<2;i++) ivec[2*ib+i]=veci[2*tx+i];\
}

#define M2L_KERNEL(m2l_kernel,m2l_kernel_core) \
extern "C" __global__ \
void m2l_kernel(int *nvec,float *ivec,float *jvec,float *kvec,float *lvec)\
{\
  int bx=blockIdx.x;\
  int tx=threadIdx.x;\
  int i,j,k,m,n,ib,jb,ij,nij,jbase,jblok,je,nms,nmk,nks,jks,jnk;\
  int mblok=nvec[2],mp=nvec[6],nrbm=nvec[7];\
  int nb=0,nc[3],nd=0,mpdnm=(4*mp*mp*mp-mp)/3;\
  float rb=scal[0],eps=1e-6;\
  float xijc=0,yijc=0,zijc=0,rh,cnmre,cnmim,dnmre,dnmim;\
  float sr,ank,ajk,ajn,fnmm,fnpm;\
  float vecd[2];\
  __shared__ int mg[NBLOK1],ng[NBLOK1];\
  __shared__ float veci[6*NBLOK1],vecj[2*NBLOK1];\
  __shared__ float ynmre[MPMAX],ynmim[MPMAX];\
  nc[0]=0;\
  rh=xijc+yijc+zijc+eps+nb+nc[0]+nd;\
  ib=bx*NBLOK1+tx;\
  nij=nvec[bx*mblok+10];\
  for(i=0;i<NBLOK1;i++){\
    ng[i]=0;\
    mg[i]=0;\
  }\
  for(n=0;n<mp;n++){\
    for(m=0;m<=n;m++){\
      nms=n*(n+1)/2+m;\
      ng[nms]=n;\
      mg[nms]=m;\
    }\
  }\
  jblok=(MPMAX+NBLOK1-1)/NBLOK1;\
  for(j=0;j<jblok-1;j++){\
    jb=j*NBLOK1+tx;\
    ynmre[jb]=kvec[2*jb+0];\
    ynmim[jb]=kvec[2*jb+2];\
    __syncthreads();\
  }\
  if(j*NBLOK1+tx<MPMAX){\
    jb=j*NBLOK1+tx;\
    ynmre[jb]=kvec[2*jb+0];\
    ynmim[jb]=kvec[2*jb+1];\
  }\
  __syncthreads();\
  for(i=0;i<6;i++) veci[6*tx+i]=0;\
  __syncthreads();\
  for(ij=0;ij<nij;ij++){\
    jbase=nvec[bx*mblok+2*ij+11];\
    je=nvec[bx*mblok+2*ij+12];\
    jb=jbase+tx;\
    for(i=0;i<2;i++) vecj[2*tx+i]=jvec[2*jb+i];\
    __syncthreads();\
    m2l_kernel_core;\
    for(i=0;i<2;i++) veci[6*tx+3*i]+=vecd[i];\
    __syncthreads();\
  }\
  for(i=0;i<2;i++) ivec[2*ib+i]=veci[6*tx+3*i];\
}

#define L2P_KERNEL(l2p_kernel,l2p_kernel_core) \
extern "C" __global__ \
void l2p_kernel(int *nvec,float *ivec,float *jvec,float *kvec,float *lvec)\
{\
  int bx=blockIdx.x;\
  int tx=threadIdx.x;\
  int i,j,k,m,n,ib,jb,jbase,jblok,nm,nms;\
  int mblok=nvec[2],mp=nvec[6],nb,nc[3],nd;\
  float rb=scal[0],xmin=scal[1],ymin=scal[2],zmin=scal[3];\
  float xic,yic,zic,xiic,yiic,ziic,r,th,ph;\
  float xx,yy=0,s2,p,pn,p1,p2,fact,ere,eim,eps=1e-6;\
  float rsre=0,rsim=0,rrre=0,rrim=0,rthre=0,rthim=0,rphre=0,rphim=0;\
  float g=0,gr=0,gth=0,gph=0,gx=0,gy=0,gz=0;\
  float bnm,bth=0;\
  __shared__ float veci[6*NBLOK1],vecj[2*NBLOK1],veck[4*NBLOK1],vecl[MPMAX];\
  r=yy+bth+rsre+rsim+rrre+rrim+rthre+rthim+rphre+rphim+g+gr+gth+gph+gx+gy+gz;\
  ib=bx*NBLOK1+tx;\
  jblok=(MPMAX+NBLOK1-1)/NBLOK1;\
  for(j=0;j<jblok-1;j++){\
    jb=j*NBLOK1+tx;\
    vecl[jb]=lvec[jb];\
    vecl[jb]=lvec[jb];\
    __syncthreads();\
  }\
  if(j*NBLOK1+tx<MPMAX){\
    jb=j*NBLOK1+tx;\
    vecl[jb]=lvec[jb];\
    vecl[jb]=lvec[jb];\
  }\
  __syncthreads();\
  nb=nvec[bx*mblok+10];\
  jbase=nvec[bx*mblok+11];\
  for(i=0;i<3;i++) nc[i]=0;\
  k=0;\
  i=1;\
  while(nb!=0){\
    j=2-k;\
    nc[j]=nc[j]+nb%2*i;\
    nb=nb/2;\
    j=k+1;\
    k=j%3;\
    if(k==0) i=i*2;\
  }\
  nd=nc[0];\
  nc[0]=nc[1];\
  nc[1]=nc[2];\
  nc[2]=nd;\
  xic=xmin+(nc[0]+0.5)*rb;\
  yic=ymin+(nc[1]+0.5)*rb;\
  zic=zmin+(nc[2]+0.5)*rb;\
  jb=jbase+tx;\
  for(i=0;i<6;i++) veci[6*tx+i]=ivec[6*ib+i];\
  for(i=0;i<2;i++) vecj[2*tx+i]=jvec[2*jb+i];\
  for(i=0;i<4;i++) veck[4*tx+i]=0.0f;\
  __syncthreads();\
  l2p_kernel_core;\
  for(i=0;i<4;i++) kvec[4*ib+i]=veck[4*tx+i];\
}

P2P_KERNEL(G_p2p_kernel,G_P2P_KERNEL_CORE);
P2P_KERNEL(Gni_p2p_kernel,GNI_P2P_KERNEL_CORE);
P2P_KERNEL(Gnj_p2p_kernel,GNJ_P2P_KERNEL_CORE);
P2M_KERNEL(G_p2m_kernel,G_P2M_KERNEL_CORE);
P2M_KERNEL(Gn_p2m_kernel,GN_P2M_KERNEL_CORE);
M2L_KERNEL(m2m_kernel,M2M_KERNEL_CORE);
M2L_KERNEL(m2l_kernel,M2L_KERNEL_CORE);
M2L_KERNEL(l2l_kernel,L2L_KERNEL_CORE);
L2P_KERNEL(G_l2p_kernel,G_L2P_KERNEL_CORE);
L2P_KERNEL(Gn_l2p_kernel,GN_L2P_KERNEL_CORE);
