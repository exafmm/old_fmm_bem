/*
number of mathematical operations (only floating point)
      operation  flo/o total
      +-* :190      1   190
      /   : 12      4    48
      sqrt:  2      4     8
      sin :  1      8     8
      cos :  3      8    24
      pow :  1     13    13
      sum               291
*/
#define G_P2M_KERNEL_CORE \
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
      s2=sqrtf((1-xx)*(1+xx));\
      fact=1;\
      pn=1;\
      rhm=1;\
      for(m=0;m<=mg[tx];m++){\
        p=pn;\
        nm=m*m+2*m;\
        bnm[nm]=rhm*veck[nm]*p;\
        p1=p;\
        p=xx*(2*m+1)*p;\
        rhm*=rh;\
        rhn=rhm;\
        for(n=m+1;n<=ng[tx];n++){\
          nm=n*n+n+m;\
          bnm[nm]=rhn*veck[nm]*p;\
          p2=p1;\
          p1=p;\
          p=(xx*(2*n+1)*p1-(n+m)*p2)/(n-m+1);\
          rhn*=rh;\
        }\
        pn=-pn*fact*s2;\
        fact=fact+2;\
      }\
      n=ng[tx];\
      m=mg[tx];\
      nm=n*n+n+m;\
      ere=__cosf(-m*be);\
      eim=__sinf(-m*be);\
      veci[2*tx+0]+=vecj[jj7+6]*bnm[nm]*ere;\
      veci[2*tx+1]+=vecj[jj7+6]*bnm[nm]*eim;\
      jj7+=7;
