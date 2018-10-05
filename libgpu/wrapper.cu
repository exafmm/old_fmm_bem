#define P2P_WRAPPER(p2p_wrapper,p2p_kernel) \
extern "C" \
void p2p_wrapper()\
{\
  dim3 block(NBLOK0);\
  dim3 grid(iblok);\
  p2p_kernel<<< grid, block >>>(nveg[idev],iveg[idev],jveg[idev],kveg[idev]);\
  CUT_CHECK_ERROR("Kernel execution failed");\
}

#define P2M_WRAPPER(p2m_wrapper,p2m_kernel) \
extern "C" \
void p2m_wrapper()\
{\
  dim3 block(NBLOK1);\
  dim3 grid(iblok);\
  p2m_kernel<<< grid, block >>>(nveg[idev],iveg[idev],jveg[idev],kveg[idev]);\
  CUT_CHECK_ERROR("Kernel execution failed");\
}

#define M2L_WRAPPER(m2l_wrapper,m2l_kernel) \
extern "C" \
void m2l_wrapper()\
{\
  dim3 block(NBLOK1);\
  dim3 grid(iblok);\
  m2l_kernel<<< grid, block >>>(nveg[idev],iveg[idev],jveg[idev],kveg[idev],lveg[idev]);\
  CUT_CHECK_ERROR("Kernel execution failed");\
}

#define L2P_WRAPPER(l2p_wrapper,l2p_kernel) \
extern "C" \
void l2p_wrapper()\
{\
  dim3 block(NBLOK1);\
  dim3 grid(iblok);\
  l2p_kernel<<< grid, block >>>(nveg[idev],iveg[idev],jveg[idev],kveg[idev],lveg[idev]);\
  CUT_CHECK_ERROR("Kernel execution failed");\
}

P2P_WRAPPER(G_p2p_gpu,G_p2p_kernel);
P2P_WRAPPER(Gni_p2p_gpu,Gni_p2p_kernel);
P2P_WRAPPER(Gnj_p2p_gpu,Gnj_p2p_kernel);
P2M_WRAPPER(G_p2m_gpu,G_p2m_kernel);
P2M_WRAPPER(Gn_p2m_gpu,Gn_p2m_kernel);
M2L_WRAPPER(m2m_gpu,m2m_kernel);
M2L_WRAPPER(m2l_gpu,m2l_kernel);
M2L_WRAPPER(l2l_gpu,l2l_kernel);
L2P_WRAPPER(G_l2p_gpu,G_l2p_kernel);
L2P_WRAPPER(Gn_l2p_gpu,Gn_l2p_kernel);
