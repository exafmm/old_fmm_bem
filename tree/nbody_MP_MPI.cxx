// Copyright(C) 2010 by 
// Tsuyoshi Hamada <hamada@progrape.jp>
// Keigo Nitadori <nitadori@margaux.astron.s.u-tokyo.ac.jp>
// Rio Yokota <rio.yokota@bristol.ac.uk>

#include <mpi.h>
using namespace std;
#define real double
#include "BHtree.h"
#include "nbody.h"
#include "delegate_alltoall.h"

static MPI_Status status;
static int local_proc_id;
static int total_proc_count;
static char MP_hostname[NMAXPROC][MPI_MAX_PROCESSOR_NAME];
static char *myname;
static MPI_Datatype MPI_FLOAT4   = 0;
static MPI_Datatype MPI_PARTICLE = 0;

int MP_myprocid()
{
  return local_proc_id;
}

int MP_proccount()
{
  return total_proc_count;
}

void MP_initialize(int *argc, char ***argv)
{
  int namelen;
  int myid, numprocs;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int prob;
  MPI_Init_thread(argc, argv, MPI_THREAD_FUNNELED, &prob);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  MPI_Get_processor_name(processor_name,&namelen);
#ifdef PRINT
  cerr << "Initialize:Myid = " << myid
       << " Myname = " << processor_name
       << " Nprocs = " << numprocs <<endl;
#endif
  total_proc_count = numprocs;
  local_proc_id = myid;
  MPI_Get_processor_name(MP_hostname[myid], &namelen);
  MPI_Allgather(MP_hostname[myid], MPI_MAX_PROCESSOR_NAME, MPI_CHAR,
          MP_hostname,       MPI_MAX_PROCESSOR_NAME, MPI_CHAR, MPI_COMM_WORLD);
  myname = MP_hostname[myid];
#ifdef PRINT
  if(MP_root()){
    for(int p=0; p<numprocs; p++){
      printf("%04d\t%s\n", p, MP_hostname[p]);
    }
  }
  for(int p=0; p<numprocs; p++){
    assert(MP_hostname[p][0] | 'g' | 'b' == 'g' | 'b');
  }
  printf("Hello, I am %s, my id is %d\n", MP_get_hostname(), MP_myprocid());
#endif
  mympi::initialize();
  MPI_Type_contiguous(4, MPI_FLOAT, &MPI_FLOAT4);
  MPI_Type_commit(&MPI_FLOAT4);
  int ss = sizeof(nbody_particle) / sizeof(double);
  assert(0 == sizeof(nbody_particle) % sizeof(double));
  MPI_Type_contiguous(ss, MPI_DOUBLE, &MPI_PARTICLE);
  MPI_Type_commit(&MPI_PARTICLE);
  MPI_Barrier(MPI_COMM_WORLD);
}

void MP_gather_sample_coords(std::vector<float4> &sample_array)
{
  static std::vector<float4> recvbuf;
  int nprocs = MP_proccount();
  static int count[NMAXPROC]; 
  static int displs[NMAXPROC+1] = {0, };
  int nsample = sample_array.size();
  MPI_Gather(&nsample, 1, MPI_INT, count, 1, MPI_INT, 0, MPI_COMM_WORLD);
  for(int i=0; i<nprocs; i++){
    displs[i+1] = displs[i] + count[i];
  }
  if(MP_root()){
    recvbuf.resize(displs[nprocs]);
  }
  MPI_Gatherv(&sample_array[0], nsample,       MPI_FLOAT4, 
              &recvbuf     [0], count, displs, MPI_FLOAT4, 0, MPI_COMM_WORLD);
  if(MP_root()){
    sample_array.swap(recvbuf);
  }
}

void MP_int_bcast(int& i)
{
  MPI_Bcast(&i,1,MPI_INT,0,MPI_COMM_WORLD);
}

void MP_double_bcast(double* data, int nwords)
{
  MPI_Bcast(data,nwords,MPI_DOUBLE,0,MPI_COMM_WORLD);
}

void MP_int_sum(int& i)
{
  int tmp;
  MPI_Reduce(&i,&tmp,1, MPI_INT, MPI_SUM,0,MPI_COMM_WORLD);
  if(local_proc_id == 0) i = tmp;
}

void MP_sum(double& r)
{
  double tmp;
  MPI_Reduce(&r,&tmp,1, MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD);
  if(local_proc_id == 0) r = tmp;
}

double MP_doublemax(double localval)
{
  double globalval = localval;
  MPI_Allreduce(&localval, &globalval,1, MPI_DOUBLE, MPI_MAX,MPI_COMM_WORLD);
  return globalval;
}

void MP_sync()
{
  MPI_Barrier(MPI_COMM_WORLD);
}

void MP_end()
{
  MPI_Finalize();
}

void exchange_local_essential_trees_alltoall_detune(nbody_particle ptcl[],
						    bhparticle bp[],
						    int nbody,
						    int nbody_max,
						    real theta2,
						    vector3 xmin[],
						    vector3 xmax[],
						    bhnode *root,
						    int &ntot){
						    int nprocs = MP_proccount();
						    int myid   = MP_myprocid();
						    static std::ofstream nullstr("/dev/null");

  static int initcall = true;
  static vector<float4> sendbuf[NMAXPROC];
  static vector<float4> recvbuf[NMAXPROC];
  if(initcall){
    initcall = false;
    try{
      for(int p=0; p<nprocs; p++){
        sendbuf[p].reserve(256);
      }
    }catch(bad_alloc){
      cout << myid << " : bad alloc 1" << endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }
  {
    RAII_Timer timer("make LET list, ", NULL, MP_root() ? cout : nullstr);
#pragma omp parallel for 
    for(int p=0; p<nprocs; p++){
      if(p != myid){
        vector3 center = 0.5 * (xmax[p] + xmin[p]);
        vector3 hlen   = 0.5 * (xmax[p] - xmin[p]);
        sendbuf[p].clear();
        root->create_essential_tree_vector(center, hlen, sendbuf[p], theta2);
      }
    }
    MP_sync();
  }
  double dtime;
  {
    RAII_Timer timer("alltoall delegate, ", &dtime, MP_root() ? cout : nullstr);
    MP_alltoallv(sendbuf, recvbuf);
  }
  int nsend_tot = 0;
  int nrecv_tot = 0;
  for(int p=0; p<nprocs; p++){
    nsend_tot += sendbuf[p].size();
    nrecv_tot += recvbuf[p].size();
  }
  assert(nrecv_tot <= nbody_max - nbody);
  MP_int_sum(nsend_tot);
  MP_int_sum(nrecv_tot);
  if(myid == 0){
       assert(nsend_tot == nrecv_tot); 
  }
#ifdef PRINT
  double bw = 2.0 * double(sizeof(float4) * nsend_tot) / dtime * 1.e-9;
  if(myid == 0) cout << "Exchanged treenodes = " 
                   << nsend_tot << endl;
  if(myid == 0) cout << "Global Bandwidth " << bw << " GB/s" << endl;
#endif

  int off = 0;
  for(int p=0; p<nprocs; p++){
    int rsize = recvbuf[p].size();
    float4 *ptr = &recvbuf[p][0];
    for(int i=0; i<rsize; i++){
      bp[nbody+off+i].set_extp(ptr[i]);
    }
    off += rsize;
  }
  ntot = nbody + off;
}

void MP_Abort(int err){
  MPI_Abort(MPI_COMM_WORLD, err);
}

bool MP_root(){
  return (0 == MP_myprocid());
}
const char *MP_get_hostname(int id){
  return MP_hostname[id];
}
const char *MP_get_hostname(){
  return MP_hostname[local_proc_id];
}

template <typename T>
void MP_alltoallv(std::vector<T> sendbuf[], std::vector<T> recvbuf[]){
  mympi::delegate_alltoallv(sendbuf, recvbuf);
}
template void MP_alltoallv<float4>(std::vector<float4> [], std::vector<float4> []);
template void MP_alltoallv<nbody_particle>(std::vector<nbody_particle> [], std::vector<nbody_particle> []);
