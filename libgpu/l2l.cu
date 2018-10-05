/*
number of mathematical operations (only floating point)
      operation  flo/o total
      +-* : 27      1    27
      /   :  2      4     8
      pow :  1     13    13
      sum                48
*/
#define L2L_KERNEL_CORE \
    rh=rb*sqrtf(3.0f)/2;\
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
    for(n=j;n<mp;n++){\
      nks=n*(n+1)/2+k;\
      jnk=(n-j)*(n-j)+n-j;\
      fnmm=1.0;\
      for(i=0;i<n-k;i++) fnmm=fnmm*(i+1);\
      fnpm=1.0;\
      for(i=0;i<n+k;i++) fnpm=fnpm*(i+1);\
      ank=pow(-1.0,n)*rsqrtf(fnmm*fnpm);\
      fnpm=1.0;\
      for(i=0;i<n-j;i++) fnpm=fnpm*(i+1);\
      ajn=pow(-1.0,n-j)/fnpm;\
      sr=ajn*ajk/ank;\
      cnmre=sr*ynmre[jnk]*pow(rh,n-j);\
      cnmim=sr*ynmim[jnk]*pow(rh,n-j);\
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
