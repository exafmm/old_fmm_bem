// Copyright(C) 2010 by 
// Tsuyoshi Hamada <hamada@progrape.jp>
// Keigo Nitadori <nitadori@margaux.astron.s.u-tokyo.ac.jp>
// Rio Yokota <rio.yokota@bristol.ac.uk>

#ifndef _NBODY_H
#define _NBODY_H

#include "nbody_particle.h"

#define NMAXSAMPLE (500 * 1000)
#define BHLISTMAX (500000)
#define LETINDEXBASE (1000000000)
#define NSIZEFACTOR     2.4
#define NSIZEFACTOR_BHP 3.0
#define NSIZEBIAS  100000

inline int arraysize_bhp(const int n)
{
  return (int) (n*NSIZEFACTOR_BHP+NSIZEBIAS);
}
void wall_init();
void MP_initialize(int *argc, char ***argv);
void MP_end();
void MP_sync();
int MP_myprocid();
int MP_proccount();
void MP_int_bcast(int&i);
void MP_int_sum(int&i);
void MP_sum(double& r);
void MP_double_bcast(double*i, int nwords);
double MP_doublemax(double localval);
void MP_Abort(int);
template <typename T>
void MP_alltoallv(std::vector<T> sendbuf[], std::vector<T> recvbuf[]);

void exchange_local_essential_trees_alltoall_detune(
  nbody_particle ptcl[],
  bhparticle bp[],
  int nbody,
  int nbody_max,
  real theta2,
  vector3 xlow[],
  vector3 xhigh[],
  bhnode *root,
  int &ntot);

bool MP_root();
const char *MP_get_hostname(int);
const char *MP_get_hostname();

// distribute.C

int set_nmaxsample(int nmaxsample);

void initialize_division(int nbody);

void setup_division(nbody_particle * pb,
		    int nbody,
		    int npdim[3],
		    vector3* xlow,
		    vector3* xhigh);
void exchange_particles_alltoall_vector(
		nbody_particle pb[],
		int &nbody,
		const int nbmax,
		const vector3 xlow[],
		const vector3 xhigh[]);
void exchange_particles_alltoall_vector_init(
		nbody_particle pb_src[],
		nbody_particle pb_dst[],
		int &nbody,
		const int nbmax,
		const vector3 xlow[],
		const vector3 xhigh[]);
void exchange_local_essential_trees(nbody_particle * pb,
				    int nbody,
				    int nbmax,
				    real theta2,
				    vector3 * xlow,
				    vector3 * xhigh,
				    bhnode * bp,
				    int & ntot);

// BHtree.h

bhnode * get_bhroot();
inline void log_allocate(const char  *message, size_t size){}

//time.C

real wall_time();
real cpu_time();

#endif





