/*
number of mathematical operations (only floating point)
      operation  flo/o total
      +-* : 23      1    23
      /   :  6      4    24
      sqrt:  1      4     4
      pow :  1     13    13
      sum                64
*/
#define GNI_P2P_KERNEL_CORE \
      dxij=veci[0]-vecj[jj7  ];\
      dyij=veci[1]-vecj[jj7+1];\
      dzij=veci[2]-vecj[jj7+2];\
      rij=rsqrtf(dxij*dxij+dyij*dyij+dzij*dzij+eps);\
      rsij=pi14*vecj[jj7+6]*rij*rij*rij;\
      veck[0]+=(dxij*veci[3]+dyij*veci[4]+dzij*veci[5])*rsij;\
      jj7+=7;
