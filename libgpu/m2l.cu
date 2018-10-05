/*
number of mathematical operations (only floating point)
      operation  flo/o total
      +-* : 27      1    27
      /   :  2      4     8
      pow :  1     13    13
      sum                48
*/
#define M2L_KERNEL_CORE \
    for(i=0;i<3;i++) nc[i]=0;\
    nb=je-1;\
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
    xijc=(nc[0]-3)*rb;\
    yijc=(nc[1]-3)*rb;\
    zijc=(nc[2]-3)*rb;\
    rh=sqrt(xijc*xijc+yijc*yijc+zijc*zijc)+eps;\
    jbase=(je-1)*mpdnm;\
    n=ng[tx];\
    m=mg[tx];\
    nms=n*(n+1)/2+m;\
    for(i=0;i<2;i++) vecd[i]=0;\
    for(k=-n;k<0;k++){\
      nks=n*(n+1)/2-k;\
      nmk=jbase+(4*n*n*n+6*n*n+5*n)/3+m*(2*n+1)+k;\
      dnmre=lvec[2*nmk+0];\
      dnmim=lvec[2*nmk+1];\
      vecd[0]+=dnmre*vecj[2*nks+0];\
      vecd[0]+=dnmim*vecj[2*nks+1];\
      vecd[1]-=dnmre*vecj[2*nks+1];\
      vecd[1]+=dnmim*vecj[2*nks+0];\
    }\
    for(k=0;k<=n;k++){\
      nks=n*(n+1)/2+k;\
      nmk=jbase+(4*n*n*n+6*n*n+5*n)/3+m*(2*n+1)+k;\
      dnmre=lvec[2*nmk+0];\
      dnmim=lvec[2*nmk+1];\
      vecd[0]+=dnmre*vecj[2*nks+0];\
      vecd[0]-=dnmim*vecj[2*nks+1];\
      vecd[1]+=dnmre*vecj[2*nks+1];\
      vecd[1]+=dnmim*vecj[2*nks+0];\
    }\
    __syncthreads();\
    for(i=0;i<2;i++) vecj[2*nms+i]=vecd[i];\
    __syncthreads();\
    j=ng[tx];\
    k=mg[tx];\
    jks=j*(j+1)/2+k;\
    for(i=0;i<2;i++) vecd[i]=0;\
    fnmm=1.0;\
    for(i=0;i<j-k;i++) fnmm=fnmm*(i+1);\
    fnpm=1.0;\
    for(i=0;i<j+k;i++) fnpm=fnpm*(i+1);\
    ajk=pow(-1.0,j)*rsqrtf(fnmm*fnpm);\
    for(n=abs(k);n<mp;n++){\
      nks=n*(n+1)/2+k;\
      jnk=(j+n)*(j+n)+j+n;\
      fnmm=1.0;\
      for(i=0;i<n-k;i++) fnmm=fnmm*(i+1);\
      fnpm=1.0;\
      for(i=0;i<n+k;i++) fnpm=fnpm*(i+1);\
      ank=pow(-1.0,n)*rsqrtf(fnmm*fnpm);\
      fnpm=1.0;\
      for(i=0;i<j+n;i++) fnpm=fnpm*(i+1);\
      ajn=pow(-1.0,j+n)/fnpm;\
      sr=pow(-1.0,j+k)*ank*ajk/ajn;\
      cnmre=sr*ynmre[jnk]/pow(rh,j+n+1);\
      cnmim=sr*ynmim[jnk]/pow(rh,j+n+1);\
      vecd[0]+=vecj[2*nks+0]*cnmre;\
      vecd[0]-=vecj[2*nks+1]*cnmim;\
      vecd[1]+=vecj[2*nks+0]*cnmim;\
      vecd[1]+=vecj[2*nks+1]*cnmre;\
    }\
    __syncthreads();\
    for(i=0;i<2;i++) vecj[2*jks+i]=vecd[i];\
    __syncthreads();\
    jbase=(je+nrbm-1)*mpdnm;\
    n=ng[tx];\
    m=mg[tx];\
    nms=n*(n+1)/2+m;\
    for(i=0;i<2;i++) vecd[i]=0;\
    for(k=-n;k<0;k++){\
      nks=n*(n+1)/2-k;\
      nmk=jbase+(4*n*n*n+6*n*n+5*n)/3+m*(2*n+1)+k;\
      dnmre=lvec[2*nmk+0];\
      dnmim=lvec[2*nmk+1];\
      vecd[0]+=dnmre*vecj[2*nks+0];\
      vecd[0]+=dnmim*vecj[2*nks+1];\
      vecd[1]-=dnmre*vecj[2*nks+1];\
      vecd[1]+=dnmim*vecj[2*nks+0];\
    }\
    for(k=0;k<=n;k++){\
      nks=n*(n+1)/2+k;\
      nmk=jbase+(4*n*n*n+6*n*n+5*n)/3+m*(2*n+1)+k;\
      dnmre=lvec[2*nmk+0];\
      dnmim=lvec[2*nmk+1];\
      vecd[0]+=dnmre*vecj[2*nks+0];\
      vecd[0]-=dnmim*vecj[2*nks+1];\
      vecd[1]+=dnmre*vecj[2*nks+1];\
      vecd[1]+=dnmim*vecj[2*nks+0];\
    }
