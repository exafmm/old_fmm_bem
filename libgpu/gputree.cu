// Copyright(C) 2010 by 
// Tsuyoshi Hamada <hamada@progrape.jp>
// Keigo Nitadori <nitadori@margaux.astron.s.u-tokyo.ac.jp>
// Rio Yokota <rio.yokota@bristol.ac.uk>

#include "vforce.h"

namespace libcunbody{
using namespace std;

__device__ float4 dev_apot(float4 xi, float4 xj, float4 apot)
{
  float dx = xj.x - xi.x;
  float dy = xj.y - xi.y;
  float dz = xj.z - xi.z;
  float eps2 = xi.w;
  float mj   = xj.w;
  float r2   = ((eps2 + dx*dx) + dy*dy) + dz*dz;
  float r1i  = rsqrt(r2);
  float r2i  = r1i*r1i;
  float mr1i = mj * r1i;
  float mr3i = mr1i * r2i;
  apot.x += dx * mr3i;
  apot.y += dy * mr3i;
  apot.z += dz * mr3i;
  apot.w -= mr1i;
  return (apot);
}

__global__ void cunbody_kernel_tree_015(
  float4 *xilist, 
  float4 *xjlist, 
  float4 *apotlist, 
  int2 *off)
{
  const int NJBLOCK = 128;
  int gid = threadIdx.x + blockDim.x * blockIdx.x;
  float4 xi = xilist[gid];
  float4 apot = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
    
  int jstart = off[gid].x;
  int jend   = off[gid].y;

  for(int jbase=jstart; jbase<jend; jbase+=NJBLOCK){
#pragma unroll 128
    for(int j=0; j<NJBLOCK; j++){
      float4 xj = xjlist[jbase + j];
      apot = dev_apot (xi, xj, apot);
    }
  }
  apotlist[gid] = apot;
}

}; // namespace libcunbody __END__
