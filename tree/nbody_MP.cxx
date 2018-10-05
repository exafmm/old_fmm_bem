// Copyright(C) 2010 by 
// Tsuyoshi Hamada <hamada@progrape.jp>
// Keigo Nitadori <nitadori@margaux.astron.s.u-tokyo.ac.jp>
// Rio Yokota <rio.yokota@bristol.ac.uk>

using namespace std;

#define real double
#include  "BHtree.h"
#include  "nbody.h"

static int local_proc_id = 0;
static int total_proc_count  = 1;

int MP_myprocid()
{
  return local_proc_id;
}

int MP_proccount()
{
  return total_proc_count;
}

void MP_initialize(int *argc,char ***argv)
{
  return;
}

double MP_doublemax(double localval)
{
  return localval;
}

void MP_sync()
{
}

void MP_end()
{
}

void MP_gather_sample_coords(std::vector<float4> &sample_array)
{
}

void MP_int_bcast(int& i)
{
}

void MP_int_sum(int& i)
{
}

void MP_sum(double& r) {}
void MP_double_bcast(double*i, int nwords) {}

void exchange_local_essential_trees_alltoall_detune(nbody_particle ptcl[],
						    bhparticle bp[],
						    int nbody,
						    int nbody_max,
						    real theta2,
						    vector3 xmin[],
						    vector3 xmax[],
						    bhnode *root,
						    int &ntot){
  ntot = nbody;
}

void MP_Abort(int err){
  exit(err);
}

bool MP_root(){
  return true;
}

const char *MP_get_hostname(int){
  return "serial";
}
const char *MP_get_hostname(){
  return "serial";
}

template <typename T>
void MP_alltoallv(std::vector<T> sendbuf[], std::vector<T> recvbuf[]){
  return;
}

template void MP_alltoallv<float4>(std::vector<float4> [], std::vector<float4> []);
template void MP_alltoallv<nbody_particle>(std::vector<nbody_particle> [], std::vector<nbody_particle> []);

namespace mympi{
  int local_size(){
    return 1;
  }
  int local_id(){
    return 0;
  }
  void sync_local(){
  }
}
