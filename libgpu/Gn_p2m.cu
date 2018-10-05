/*
number of mathematical operations (only floating point)
      operation  flo/o total
      +-* :326      1   326
      /   : 40      4   160
      sqrt:  2      4     8
      sin :  3      8    24
      cos :  3      8    24
      pow :  2     13    26
      sum               568
*/
#define GN_P2M_KERNEL_CORE \
      float yy;\
      float brerh,breal,brebe,bimrh,bimal,bimbe;\
      float brexd,breyd,brezd,bimxd,bimyd,bimzd;\
      __shared__ float bth[MPMAX];\
      xjjc=vecj[jj7+0]-xjc;\
      yjjc=vecj[jj7+1]-yjc;\
      zjjc=vecj[jj7+2]-zjc;\
      rh=sqrtf(xjjc*xjjc+yjjc*yjjc+zjjc*zjjc)+eps;\
      al=acosf(zjjc/rh);\
      if(abs(xjjc)+abs(yjjc)<eps){\
        be=0;\
      }\
      else if(abs(xjjc)<eps){\
        be=yjjc/abs(yjjc)*M_PI*0.5;\
      }\
      else if(xjjc>0){\
        be=atanf(yjjc/xjjc);\
      }\
      else{\
        be=atanf(yjjc/xjjc)+M_PI;\
      }\
      xx=__cosf(al);\
      yy=__sinf(al);\
      if(fabs(yy)<eps) yy=1/eps;\
      s2=sqrtf((1-xx)*(1+xx));\
      fact=1;\
      pn=1;\
      for(m=0;m<=mg[tx];m++){\
        p=pn;\
        nm=m*m+2*m;\
        bnm[nm]=veck[nm]*p;\
        p1=p;\
        p=xx*(2*m+1)*p;\
        bth[nm]=veck[nm]*(p-(m+1)*xx*p1)/yy;\
        for(n=m+1;n<=ng[tx];n++){\
          nm=n*n+n+m;\
          bnm[nm]=veck[nm]*p;\
          p2=p1;\
          p1=p;\
          p=(xx*(2*n+1)*p1-(n+m)*p2)/(n-m+1);\
          bth[nm]=veck[nm]*((n-m+1)*p-(n+1)*xx*p1)/yy;\
        }\
        pn=-pn*fact*s2;\
        fact=fact+2;\
      }\
      n=ng[tx];\
      m=mg[tx];\
      nms=n*(n+1)/2+m;\
      nm=n*n+n+m;\
      ere=__cosf(-m*be);\
      eim=__sinf(-m*be);\
      rhm=powf(rh,n-1);\
      rhn=rhm*rh;\
      brerh=n*rhm*bnm[nm]*ere;\
      breal=rhn*bth[nm]*ere;\
      brebe=m*rhn*bnm[nm]*eim;\
      bimrh=n*rhm*bnm[nm]*eim;\
      bimal=rhn*bth[nm]*eim;\
      bimbe=-m*rhn*bnm[nm]*ere;\
      brexd=__sinf(al)*__cosf(be)*brerh+__cosf(al)*__cosf(be)/rh*breal-\
      __sinf(be)/rh/yy*brebe;\
      breyd=__sinf(al)*__sinf(be)*brerh+__cosf(al)*__sinf(be)/rh*breal+\
      __cosf(be)/rh/yy*brebe;\
      brezd=__cosf(al)*brerh-__sinf(al)/rh*breal;\
      bimxd=__sinf(al)*__cosf(be)*bimrh+__cosf(al)*__cosf(be)/rh*bimal-\
      __sinf(be)/rh/yy*bimbe;\
      bimyd=__sinf(al)*__sinf(be)*bimrh+__cosf(al)*__sinf(be)/rh*bimal+\
      __cosf(be)/rh/yy*bimbe;\
      bimzd=__cosf(al)*bimrh-__sinf(al)/rh*bimal;\
      veci[2*tx+0]+=vecj[jj7+6]*vecj[jj7+3]*brexd;\
      veci[2*tx+0]+=vecj[jj7+6]*vecj[jj7+4]*breyd;\
      veci[2*tx+0]+=vecj[jj7+6]*vecj[jj7+5]*brezd;\
      veci[2*tx+1]+=vecj[jj7+6]*vecj[jj7+3]*bimxd;\
      veci[2*tx+1]+=vecj[jj7+6]*vecj[jj7+4]*bimyd;\
      veci[2*tx+1]+=vecj[jj7+6]*vecj[jj7+5]*bimzd;\
      jj7+=7;
