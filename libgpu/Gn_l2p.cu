/*
number of mathematical operations (only floating point)
      operation  flo/o total
      +-* : 37      1    37
      /   :  3      4    12
      sin :  1      8     8
      cos :  1      8     8
      pow :  8     13   104
      sum               169
*/
#define GN_L2P_KERNEL_CORE \
      xiic=veci[6*tx+0]-xic;\
      yiic=veci[6*tx+1]-yic;\
      ziic=veci[6*tx+2]-zic;\
      r=sqrtf(xiic*xiic+yiic*yiic+ziic*ziic)+eps;\
      th=acosf(ziic/r);\
      if(abs(xiic)+abs(yiic)<eps){\
        ph=0;\
      }\
      else if(abs(xiic)<eps){\
        ph=yiic/abs(yiic)*M_PI*0.5;\
      }\
      else if(xiic>0){\
        ph=atanf(yiic/xiic);\
      }\
      else{\
        ph=atanf(yiic/xiic)+M_PI;\
      }\
      gr=0;\
      gth=0;\
      gph=0;\
      xx=__cosf(th);\
      yy=__sinf(th);\
      if(fabs(yy)<eps) yy=1/eps;\
      s2=sqrtf((1-xx)*(1+xx));\
      fact=1;\
      pn=1;\
      for(m=0;m<mp;m++){\
        p=pn;\
        nm=m*m+2*m;\
        nms=m*(m+1)/2+m;\
        ere=__cosf(m*ph);\
        if(m==0) ere=0.5;\
        eim=__sinf(m*ph);\
        bnm=vecl[nm]*p;\
        p1=p;\
        p=xx*(2*m+1)*p;\
        bth=vecl[nm]*(p-(m+1)*xx*p1)/yy;\
        rrre=m*powf(r,m-1)*bnm*ere;\
        rthre=powf(r,m)*bth*ere;\
        rphre=-m*powf(r,m)*bnm*eim;\
        rrim=m*powf(r,m-1)*bnm*eim;\
        rthim=powf(r,m)*bth*eim;\
        rphim=m*powf(r,m)*bnm*ere;\
        gr+=2*(rrre*vecj[2*nms+0]-rrim*vecj[2*nms+1]);\
        gth+=2*(rthre*vecj[2*nms+0]-rthim*vecj[2*nms+1]);\
        gph+=2*(rphre*vecj[2*nms+0]-rphim*vecj[2*nms+1]);\
        for(n=m+1;n<mp;n++){\
          nm=n*n+n+m;\
          nms=n*(n+1)/2+m;\
          bnm=vecl[nm]*p;\
          p2=p1;\
          p1=p;\
          p=(xx*(2*n+1)*p1-(n+m)*p2)/(n-m+1);\
          bth=vecl[nm]*((n-m+1)*p-(n+1)*xx*p1)/yy;\
          rrre=n*powf(r,n-1)*bnm*ere;\
          rthre=powf(r,n)*bth*ere;\
          rphre=-m*powf(r,n)*bnm*eim;\
          rrim=n*powf(r,n-1)*bnm*eim;\
          rthim=powf(r,n)*bth*eim;\
          rphim=m*powf(r,n)*bnm*ere;\
          gr+=2*(rrre*vecj[2*nms+0]-rrim*vecj[2*nms+1]);\
          gth+=2*(rthre*vecj[2*nms+0]-rthim*vecj[2*nms+1]);\
          gph+=2*(rphre*vecj[2*nms+0]-rphim*vecj[2*nms+1]);\
        }\
        pn=-pn*fact*s2;\
        fact=fact+2;\
      }\
      gx=__sinf(th)*__cosf(ph)*gr+__cosf(th)*__cosf(ph)/r*gth-\
        __sinf(ph)/r/yy*gph;\
      gy=__sinf(th)*__sinf(ph)*gr+__cosf(th)*__sinf(ph)/r*gth+\
        __cosf(ph)/r/yy*gph;\
      gz=__cosf(th)*gr-__sinf(th)/r*gth;\
      veck[tx]-=0.25/M_PI*(gx*veci[6*tx+3]+gy*veci[6*tx+4]+gz*veci[6*tx+5]);
