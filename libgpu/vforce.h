// Copyright(C) 2010 by 
// Tsuyoshi Hamada <hamada@progrape.jp>
// Keigo Nitadori <nitadori@margaux.astron.s.u-tokyo.ac.jp>
// Rio Yokota <rio.yokota@bristol.ac.uk>

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include "cuda.h"
#include "cutil-mp.h"
#include <unistd.h>

namespace libcunbody{
using namespace std;

__device__ float4 dev_apot(float4 xi, float4 xj, float4 apot);
__global__ void cunbody_kernel_tree_015(float4 *xilist, float4 *xjlist, float4 *apotlist, int2 *offlist);
__global__ void cunbody_kernel_tree_016(float4 *xilist, float4 *xjlist, float4 *apotlist, int2 *offlist);

#define KERNEL_TYPE (0x015)
#define ACC_TYPE (1)
#define MAX_WALK (256)
#define NIALLOC ( 2000 * MAX_WALK)
#define NJALLOC (20000 * MAX_WALK)
#define CUDA_MALLOC_TYPE (0)
#define CUDA_MALLOC_HOST_TYPE (0)
#define DEV_OPEN_STRATEGY (1)
#define NTHRE (128)
#define IBLOCK (1)
//#define ENABLE_DITAILED_OPENMSG

class cunbody
{
private:
  bool is_open;
  int total_proc_count;      // total number of proccess in the run
  int host_pid;      // host process id
  int host_tid;      // host thread Id
  int num_threads;
  int num_threads_per_gpu;
  float4 *dev_xi;
  float4 *dev_xj;
  float4 *dev_apot;
  int2 *dev_off, *host_off;
  unsigned int isize;
  unsigned int jsize;
  unsigned int max_isize;
  unsigned int max_jsize;

  void dev_check(){
    int ndev;
    CUDA_SAFE_CALL(cudaGetDeviceCount(&ndev));
    if(ndev == 0){
      fprintf(stdout, "ndev = %d @ %s|%d\n", ndev, __FILE__, __LINE__);
      fprintf(stdout, "There is no GPUs.\n");
      exit(-1);
    }else{
      int dev = ((host_pid * num_threads + host_tid) / num_threads_per_gpu) % ndev;
      CUDA_SAFE_CALL(cudaSetDevice(dev));
      int dev0;
      cudaGetDevice(&dev0);
      assert(dev == dev0);
#ifdef PRINT
      printf("[%03d][%01d]  : cudaSetDevice to %d-th GPU\n", host_pid, host_tid, dev+1);
#endif
      cudaDeviceProp deviceProp;
      CUDA_SAFE_CALL(cudaGetDeviceProperties(&deviceProp, dev));
      if (deviceProp.major == 9999 && deviceProp.minor == 9999){
        printf("There is no device supporting CUDA.\n");
      }
#if defined(ENABLE_DITAILED_OPENMSG)
      printf("[%d]  : Major revision number:  %d\n",        host_pid, deviceProp.major);
      printf("[%d]  : Minor revision number:  %d\n",        host_pid, deviceProp.minor);
      printf("[%d]  : core clock rate:  %.2f GHz\n",        host_pid, deviceProp.clockRate * 1e-6f);
#  if  (CUDART_VERSION >= 2000)
      printf("[%d]  : Number of cores:  %d\n",              host_pid, 8 * deviceProp.multiProcessorCount);
      printf("[%d]  : Concurrent copy and execution: %s\n", host_pid, deviceProp.deviceOverlap ? "Yes" : "No");
#  endif
#endif
    }
  }

  void dev_open(int hpid, int htid, int nthreads, int nthreads_per_gpu) {
    host_pid = hpid;
    host_tid = htid;
    num_threads = nthreads;
    num_threads_per_gpu = nthreads_per_gpu;
#ifdef PRINT
    fprintf(stdout, "[%d] CUNBODY-1 library: rev.nitadori20091029 (=^.^=)/\n", host_pid);
#endif
#if defined(ENABLE_DITAILED_OPENMSG)
    fprintf(stdout, "[%d]  open GPU by host thread %d\n", host_pid, host_tid);
    fprintf(stdout, "[%d]  : KERNEL_TYPE %03x\n",     host_pid, KERNEL_TYPE);
    fprintf(stdout, "[%d]  : MAX_WALK %d\n",          host_pid, MAX_WALK);
    fprintf(stdout, "[%d]  : NTHRE    %d\n",          host_pid, NTHRE);
    fprintf(stdout, "[%d]  : IBLOCK   %d\n",          host_pid, IBLOCK);
    fprintf(stdout, "[%d]  : ACC_TYPE %d\n",          host_pid, ACC_TYPE);
    fprintf(stdout, "[%d]  : CUDA_MALLOC_TYPE %d\n",  host_pid, CUDA_MALLOC_TYPE);
    fprintf(stdout, "[%d]  : CUDA_MALLOC_HOST_TYPE %d\n", host_pid, CUDA_MALLOC_HOST_TYPE);
    fprintf(stdout, "[%d]  : DEV_OPEN_STRATEGY %d\n", host_pid, DEV_OPEN_STRATEGY);
#endif

    this->dev_check();
    max_isize = NIALLOC;
    max_jsize = NJALLOC;

#ifdef PRINT
    printf("[%d-proc, %d-thread] ******** cudaMalloc in total at dev_open() : %d MB (j=%d, i=%d)(j=%d MB, i=%d MB)\n",
      host_pid, host_tid,
      (
        (NTHRE + max_isize) * sizeof(float4) + 
        max_isize * sizeof(float4) + 
        max_jsize * sizeof(float4) + 
        (MAX_WALK + 1) * sizeof(unsigned int) + 
        (MAX_WALK + 1) * sizeof(unsigned int)
      )/(1024*1024),
      max_jsize, max_isize,
      (max_jsize * sizeof(float4))/(1024*1024), // <----------- only j
      ((NTHRE + max_isize) * sizeof(float4) +
      max_isize * sizeof(float4))/(1024*1024)  // <----------- only i
    );
#endif

    CUDA_SAFE_CALL(cudaMalloc((void **)&dev_xi,   (NTHRE + max_isize) * sizeof(float4)));
    CUDA_SAFE_CALL(cudaMalloc((void **)&dev_apot, (NTHRE + max_isize) * sizeof(float4)));
    CUDA_SAFE_CALL(cudaMalloc((void **)&dev_off,  (NTHRE + 1 + max_isize/IBLOCK) * sizeof(int2)));
    CUDA_SAFE_CALL(cudaMalloc((void **)&dev_xj,   max_jsize * sizeof(float4)));
    CUDA_SAFE_CALL(cudaHostAlloc((void **)&host_off,  (1 + max_isize/IBLOCK) * sizeof(int2), cudaHostAllocDefault));
    is_open = true;
  }

public:
  void dev_close(void) {
    if(dev_xi) CUDA_SAFE_CALL(cudaFree(dev_xi));
    dev_xi = NULL;
    if(dev_xj) CUDA_SAFE_CALL(cudaFree(dev_xj));
    dev_xj = NULL;
    if(dev_apot) CUDA_SAFE_CALL(cudaFree(dev_apot));
    dev_apot = NULL;
    if(dev_off) CUDA_SAFE_CALL(cudaFree(dev_off));
    dev_off = NULL;
    if(host_off) CUDA_SAFE_CALL(cudaFreeHost(host_off));
  }

  cunbody() {
    is_open = false;
    total_proc_count = -1; // this is used only Wait-cudaMalloc
    host_pid = -1;
    host_tid = -1;
    max_isize = 0;
    max_jsize = 0;
    dev_xi = NULL;
    dev_xj = NULL;
    dev_apot = NULL;
    dev_off = NULL;
    host_off = NULL;
  }

  ~cunbody() {
    this->dev_close();
    dev_xi = NULL;
    dev_xj = NULL;
    dev_apot = NULL;
    dev_off = NULL;
    host_off = NULL;
    max_isize = 0;
    max_jsize = 0;
    host_tid = -1;
    host_pid = -1;
    is_open = false;
  }

  void set_total_proc_count(int n){ this->total_proc_count = n; }

  void vforce_open(int host_pid, int host_tid, int nthreads, int nthreads_per_gpu)
  {
    if(is_open == false) this->dev_open(host_pid, host_tid, nthreads, nthreads_per_gpu);
  }

  void vforce_1st(
    int host_pid, 
    int host_tid, 
    float4 xilist[], 
    float4 xjlist[], 
    unsigned int ioff[], 
    unsigned int joff[], 
    unsigned int nwalk)
  {
    isize = ioff[nwalk];
    jsize = joff[nwalk];

    if(is_open == false) this->dev_open(host_pid, host_tid, num_threads, num_threads_per_gpu);
    if(isize > max_isize){
      int isize_bak = isize;
      isize = (int)(isize*1.1);
      max_isize = isize_bak;
      CUDA_SAFE_CALL(cudaFree(dev_xi));
      CUDA_SAFE_CALL(cudaFree(dev_apot));
      CUDA_SAFE_CALL(cudaFree(dev_off));
      CUDA_SAFE_CALL(cudaFreeHost(host_off));
      CUDA_SAFE_CALL(cudaMalloc((void **)&dev_xi,   (NTHRE + isize) * sizeof(float4)));
      CUDA_SAFE_CALL(cudaMalloc((void **)&dev_apot, (NTHRE + isize) * sizeof(float4)));
      CUDA_SAFE_CALL(cudaMalloc((void **)&dev_off,  (NTHRE + 1 + isize/IBLOCK) * sizeof(int2)));
      CUDA_SAFE_CALL(cudaHostAlloc((void **)&host_off,  (1 + max_isize/IBLOCK) * sizeof(int2), cudaHostAllocDefault));
      isize = isize_bak;
    }

    if(jsize > max_jsize){
      max_jsize = jsize;
      CUDA_SAFE_CALL(cudaFree(dev_xj));
      CUDA_SAFE_CALL(cudaMalloc((void **)&dev_xj,   jsize * sizeof(float4)));
    }

    for(int iw=0; iw<nwalk; iw++){
      assert(ioff[iw]%IBLOCK == 0);
      for(int i=ioff[iw]/IBLOCK; i<ioff[iw+1]/IBLOCK; i++){
        host_off[i] = make_int2(joff[iw], joff[iw+1]);
      }
    }
    int nirun;
    for(nirun=ioff[nwalk]/IBLOCK; nirun%NTHRE;  nirun++){
      host_off[nirun] = make_int2(0, 0);
    }
    assert(nirun*IBLOCK <= max_isize);
    assert(nirun%NTHRE == 0);

    CUDA_SAFE_CALL(cudaMemcpy(dev_xi, xilist, isize * sizeof(float4), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(dev_off, host_off, nirun * sizeof(int2), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(dev_xj, xjlist, jsize * sizeof(float4), cudaMemcpyHostToDevice));

#if (KERNEL_TYPE == 0x013)
    dim3 grid(nwalk);
    dim3 threads(NTHRE);
    cunbody_kernel_tree_013 <<< grid, threads >>> (dev_xi, dev_xj, dev_apot, dev_ioff, dev_joff);
#elif (KERNEL_TYPE == 0x014)
    dim3 grid(nwalk);
    dim3 threads(NTHRE);
    cunbody_kernel_tree_014 <<< grid, threads >>> (dev_xi, dev_xj, dev_apot, dev_ioff, dev_joff);
#elif (KERNEL_TYPE == 0x015)
    dim3 grid(nirun/NTHRE);
    dim3 threads(NTHRE);
    cunbody_kernel_tree_015 <<< grid, threads >>> (dev_xi, dev_xj, dev_apot, dev_off);
#elif (KERNEL_TYPE == 0x016)
    dim3 grid(nirun/NTHRE);
    dim3 threads(NTHRE);
    cunbody_kernel_tree_016 <<< grid, threads >>> (dev_xi, dev_xj, dev_apot, dev_off);
#else
    compilation_fail
#endif
  }

  void vforce_2nd(
    int host_pid, 
    int host_tid, 
    float4 apotlist[])
  {
    CUT_CHECK_ERROR("KERNEL EXECUTION FAILED");
    CUDA_SAFE_CALL(cudaMemcpy(apotlist, dev_apot, isize * sizeof(float4), cudaMemcpyDeviceToHost));
  }

}; // class cunbody __END__

}; // namespace libcunbody __END__


#define MAX_OMP_THRE (8)
static libcunbody::cunbody cunObj[MAX_OMP_THRE];

extern "C"
{
  void vforce_1st(int proc_id, int host_tid, float4 xilist[], float4 xjlist[], unsigned int ioff[], unsigned int joff[], unsigned int nwalk){
    using namespace libcunbody;
    cunObj[host_tid].vforce_1st(proc_id, host_tid, xilist, xjlist, ioff, joff, nwalk);
  }

  void vforce_2nd(int proc_id, int host_tid, float4 apotlist[]){
    using namespace libcunbody;
    cunObj[host_tid].vforce_2nd(proc_id, host_tid, apotlist);
  }

  void vforce_mp(int proc_id, int host_tid, float4 xilist[], float4 xjlist[], float4 apotlist[], unsigned int ioff[], unsigned int joff[], unsigned int nwalk)
  {
    vforce_1st(proc_id, host_tid, xilist, xjlist, ioff, joff, nwalk);
    vforce_2nd(proc_id, host_tid, apotlist);
  }

  void vforce_open(int proc_id, int host_tid, int nthreads, int nthreads_per_gpu)
  {
    using namespace std;
    using namespace libcunbody;
    cunObj[host_tid].vforce_open(proc_id, host_tid, nthreads, nthreads_per_gpu);
  }

  void vforce_wait_randomly()
  {
  }

  void vforce_wait(double total_sec, int nproc, int pid, int tid)
  {
    int id = (pid/4) + (pid%4) * (nproc/4);
    double usec = total_sec * 1.0e+6 / (nproc*4.0);
    double nwait = usec * id;
    usleep((unsigned long) nwait);
  }

  void vforce_close(int host_tid){
#ifdef PRINT
    printf("closing cunObj[%d]\n", host_tid);
#endif
    cunObj[host_tid].dev_close();
  }

}
