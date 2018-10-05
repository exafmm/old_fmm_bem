// Copyright(C) 2010 by 
// Tsuyoshi Hamada <hamada@progrape.jp>
// Keigo Nitadori <nitadori@margaux.astron.s.u-tokyo.ac.jp>
// Rio Yokota <rio.yokota@bristol.ac.uk>

#include <mpi.h>
#include <assert.h>
#include <boost/multi_array.hpp>

#ifndef BHTREE_H
struct float4{
  typedef float  v4sf __attribute__ ((vector_size(16)));
  v4sf val;
};

struct nbody_particle{
  typedef float  v4sf __attribute__ ((vector_size(16)));
  v4sf dum[6];
};
#endif

namespace mympi{

static int nproc_pre_node = 4;

static MPI_Comm comm_world;
static MPI_Comm comm_global;
static MPI_Comm comm_local;

static int myid_world;
static bool isroot;
static int myid_global;
static int myid_local;

static int size_world;
static int size_global;
static int size_local;

const int MAX_WORLD  = 2048; // total procs
const int MAX_GLOBAL =  256; // num of node
const int MAX_LOCAL  =    8; // nprocs per node

static MPI_Datatype MPI_FLOAT3   = 0;
static MPI_Datatype MPI_FLOAT4   = 0;
static MPI_Datatype MPI_PARTICLE = 0;

static int firstcall = true;

static void initialize(const int nproc = 1, const MPI_Comm comm = MPI_COMM_WORLD)
{
  nproc_pre_node = nproc;

  //initialize communicators
  comm_world = comm;
  MPI_Comm_rank(comm_world, &myid_world);
  MPI_Comm_size(comm_world, &size_world);
  isroot = (myid_world == 0);
  assert(size_world % nproc_pre_node == 0);
  assert(size_world <= MAX_WORLD);

  MPI_Comm_split(
    comm_world, 
    myid_world / nproc_pre_node, // color
    myid_world % nproc_pre_node, // key
    &comm_local);
  MPI_Comm_rank(comm_local, &myid_local);
  MPI_Comm_size(comm_local, &size_local);
  assert(size_local <= MAX_LOCAL);

  MPI_Comm_split(
    comm_world, 
    myid_world % nproc_pre_node, // color
    myid_world / nproc_pre_node, // key
    &comm_global);

  MPI_Comm_rank(comm_global, &myid_global);
  MPI_Comm_size(comm_global, &size_global);
  assert(size_global <= MAX_GLOBAL);

  assert(size_world == size_global * size_local);

  //initialize datatypes
  MPI_Type_contiguous(3, MPI_FLOAT, &MPI_FLOAT3);
  MPI_Type_commit(&MPI_FLOAT3);

  MPI_Type_contiguous(4, MPI_FLOAT, &MPI_FLOAT4);
  MPI_Type_commit(&MPI_FLOAT4);

  int ss = sizeof(nbody_particle) / sizeof(double);
  assert(0 == sizeof(nbody_particle) % sizeof(double));
  MPI_Type_contiguous(ss, MPI_DOUBLE, &MPI_PARTICLE);
  MPI_Type_commit(&MPI_PARTICLE);

  firstcall = false;
}
template <typename T>
MPI_Datatype datatype();

template <> MPI_Datatype datatype<int> (){ return MPI_INT; }
template <> MPI_Datatype datatype<float4> (){ return MPI_FLOAT4; }
template <> MPI_Datatype datatype<nbody_particle> (){ return MPI_PARTICLE; }

#include <vector>

typedef boost::multi_array<int, 2> int2D;
template <typename T>
static void gather_vector(vector<T> dst[], vector<T> src[]){
  // NOT YET OPTIMIEZED
  static int count[MAX_LOCAL];
  static int off[MAX_LOCAL+1] = {0, };
  for(int p=0; p<size_world; p++){
    int ssize = src[p].size();
    MPI_Gather(&ssize, 1, MPI_INT, count, 1, MPI_INT, 0, comm_local);
    for(int pp=0; pp<size_local; pp++){
      off[pp+1] = off[pp] + count[pp];
    }
    if(myid_local == 0){
      dst[p].resize(off[size_local]);
    }
    MPI_Gatherv(&src[p][0], ssize, datatype<T>(), 
                &dst[p][0], count, off, datatype<T>(), 0, comm_local);
  }
}

template <typename T>
static void connect_vector(vector<T> dst[], vector<T> src[]){
  for(int ii=0, iii=0; ii<size_global; ii++){
    dst[ii].clear();
    for(int i=0; i<size_local; i++, iii++){
      dst[ii].insert(dst[ii].end(), src[iii].begin(), src[iii].end());
    }
  }
}

template <typename T>
static void alltoall_safe(vector<T> rbuf[], vector<T> sbuf[]){
  for(int dist=1; dist<size_global; dist++){
    int src = (size_global + myid_global - dist) % size_global;
    int dst = (size_global + myid_global + dist) % size_global;
    int scount = sbuf[dst].size();
    int rcount;
    MPI_Status stat;
    MPI_Sendrecv(&scount, 1, MPI_INT, dst, 0,
             &rcount, 1, MPI_INT, src, 0, comm_global, &stat);
    rbuf[src].resize(rcount);
    MPI_Sendrecv(&sbuf[dst][0], scount, datatype<T>(), dst, 1,
             &rbuf[src][0], rcount, datatype<T>(), src, 1, comm_global, &stat);
  }
  rbuf[myid_global] = sbuf[myid_global]; // just copy it
}

template <typename T>
static void scatter_vector(vector<T> dst[], vector<T> src[], int2D &counts){
  // NOT YET OPTIMIEZED
  for(int ii=0; ii<size_global; ii++){
    int *count = &counts[ii][0];
    static int off[MAX_LOCAL+1] = {0, }; 
    for(int pp=0; pp<size_local; pp++){
      off[pp+1] = off[pp] + count[pp];
    }
    int rcount = count[myid_local];
    dst[ii].resize(rcount);
    MPI_Scatterv(&src[ii][0], count, off, datatype<T>(),
             &dst[ii][0], rcount, datatype<T>(), 0, comm_local);
  }
}

template <typename T>
int delegate_alltoallv(
    vector<T> sendvec[],
    vector<T> recvvec[]){
  if(firstcall) initialize();

  // phase-1 (local comm)
  static vector<T> sendbuf1[MAX_WORLD];
  {
#ifdef PRINT
    RAII_Timer timer("phase-1 : ", NULL, std::cout, isroot);
#endif
    gather_vector(sendbuf1, sendvec);
    if(myid_local == 0){
      MPI_Barrier(comm_global);
    }
  }

  // phase-2 (global comm)
  int2D scount(boost::extents[size_global][size_local]);
  int2D rcount(boost::extents[size_global][size_local]);
  {
#ifdef PRINT
    RAII_Timer timer("phase-2 : ", NULL, std::cout, isroot);
#endif
    if(myid_local == 0){
      for(int ii=0, iii=0; ii<size_global; ii++){
        for(int i=0; i<size_local; i++, iii++){
          scount[ii][i] = sendbuf1[iii].size();
        }
      }
      MPI_Alltoall(&scount[0][0], size_local, MPI_INT,
             &rcount[0][0], size_local, MPI_INT, comm_global);
    }
    // local comm
    MPI_Bcast(&rcount[0][0], size_world, MPI_INT, 0, comm_local);
  }

  // phase-3 (no comm)
  static vector<T> sendbuf2[MAX_GLOBAL];
  {
#ifdef PRINT
    RAII_Timer timer("phase-3 : ", NULL, std::cout, isroot);
#endif
    connect_vector(sendbuf2, sendbuf1);
  }

  // phase-4 (global comm)
  static vector<T> recvbuf[MAX_GLOBAL];
  {
    RAII_Timer timer("phase-4 : ", NULL, std::cout, isroot);
    if(myid_local == 0){
#ifdef PRINT
      if(isroot) std::cout << "phase-4 : safe version" << std::endl;
#endif
      alltoall_safe(recvbuf, sendbuf2);
      MPI_Barrier(comm_global);
    }
    MPI_Barrier(comm_local);
  }

  // phase-5 (local comm)
  {
#ifdef PRINT
    RAII_Timer timer("phase-5 : ", NULL, std::cout, isroot);
#endif
    for(int i=0; i<size_global; i++){
      recvvec[i].clear();
    }
    scatter_vector(recvvec, recvbuf, rcount);
  }
  MPI_Barrier(comm_world);
  return 0;
}

int local_size(){
  return size_local;
}
int local_id(){
  return myid_local;
}
void sync_local(){
  MPI_Barrier(comm_local);
}
} // namespace mympi
