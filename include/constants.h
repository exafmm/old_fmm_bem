#define ALLOCATE_IN_MAIN
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sys/time.h>

// Adjust these parameters
const int mp      = 9;          // order of expansion in FMM
const int ngpu    = 4;           // number of GPUs per node

// The parameters below rarely require modification
const int mprint  = 0;           // 0: don't, 1: print memory usage
const int nwmax   = 50000;       // max of wall elements
const int npmax   = 1000000;     // max of particles
const int mpmax   = mp*mp;       // max of FMM coefficients
const int mpsym   = mp*(mp+1)/2; // max of FMM coefficients (symmetry)
const int mpcmp   = mp;          // max of FMM coefficients (compressed)
const int nebm    = 189;         // max of interacting boxes
const int nrbm    = 512;         // max of relative box positioning
const int nmpi    = 2;           // 1: alltoall 2: send recv 3: isend irecv
const int nsplit  = 4;           // size of local MPI communicator
const int mpim    = 500000;      // max of MPI buffer
const int nimax   = 300000;      // max of GPU i particle memory
const int njmax   = 200000;      // max of GPU j particle memory
const int nblok0  = 128;         // size of GPU thread block P2P
const int nblok1  = 64;          // size of GPU thread block B2B
const int nblok2  = 4;           // size of CPU thread block
const int ins     = 50;          // insertion threshold for quicksort
const int ngl     = 3;           // Gauss triangle 1,3,4,6,7,9,12,13
const float eps   = 1e-6;        // single precision epsilon
const float pi    = M_PI;

#ifdef MAIN
int ncheck;                  // check point switch
int nbmax;                   // total of FMM boxes @ lmax
int nbne;                    // non-empty FMM boxes @ lmax
int nbnet;                   // total of nbne for all levels
int nbned;
int lmax;                    // number of FMM box divisions
int m2lrp2p;                 // ratio of kernel weight between m2l and p2p
int nsub;                    // number of boxes in subtree
int nwn;                     // number of boundary nodes
int nwe;                     // number of bounary elements
int nch;                     // number of exterior charges
int nprocs;                  // number of processes in MPI
int myrank;                  // rank of node
int maxfour;                 // sum of Fourier modes
int maxwave;                 // sum of plane wave modes
int maxexp;                  // max of maxfour and maxwave +1
int nmax;                    // max of Fourier modes
float mem;                   // memory usage
float umem;                  // total memory usage
float xmin;                  // xmin of FMM box
float ymin;                  // ymin of FMM box
float zmin;                  // zmin of FMM box
float rd;                    // length of FMM box
double lambda;               // Helmholtz wave number
double scale;                // scaling of Bessel function
double tic,toc;              // timer

double get_time(void) {
  struct timeval tv;
  struct timezone tz;
  gettimeofday(&tv, &tz);
  return ((double)(tv.tv_sec+tv.tv_usec*1.0e-6));
}

#else
extern int ncheck;
extern int nbmax;
extern int nbne;
extern int nbnet;
extern int nbned;
extern int lmax;
extern int m2lrp2p;
extern int nsub;
extern int nwn;
extern int nwe;
extern int nch;
extern int nprocs;
extern int myrank;
extern int maxfour;
extern int maxwave;
extern int maxexp;
extern int nmax;
extern float mem;
extern float umem;
extern float xmin;
extern float ymin;
extern float zmin;
extern float rd;
extern double lambda;
extern double scale;
extern double tic,toc;

extern double get_time(void);
#endif

template<typename write_type>
void binary_write(std::fstream& fid, write_type *data, const int n) {
  int ibuf,i,j;
  const int nbuf=200000/sizeof(write_type);
  write_type buffer[nbuf];

  for( i=0; i<(n+nbuf-1)/nbuf; i++ ) {
    ibuf = std::min(nbuf,n-i*nbuf);
    for( j=0; j<ibuf; j++ ) buffer[j] = data[j+i*nbuf];
    fid.write((char *)(&buffer),sizeof(write_type)*ibuf);
  }
}

template<typename read_type>
void binary_read(std::fstream& fid, read_type *data, const int n) {
  int ibuf,i,j;
  const int nbuf=200000/sizeof(read_type);
  read_type buffer[nbuf];

  for( i=0; i<(n+nbuf-1)/nbuf; i++ ) {
    ibuf = std::min(nbuf,n-i*nbuf);
    fid.read((char *)(&buffer),sizeof(read_type)*ibuf);
    for( j=0; j<ibuf; j++ ) data[j+i*nbuf] = buffer[j];
  }
}
