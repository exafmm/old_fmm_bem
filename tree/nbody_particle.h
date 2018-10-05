#ifndef  NBODY_PARTICLE_H
#define  NBODY_PARTICLE_H
/*-----------------------------------------------------------------------------
 *  nbody-particle : basic class for simple nbody implementation
 *  J. Makino 1998/11/29
 *-----------------------------------------------------------------------------
 */

int MP_myprocid();
const char *MP_get_hostname();
void MP_Abort(int);

#define NMAXPROC 1024

#include <algorithm>
#include <cassert>
#include <fstream>
#include <malloc.h>
#include <stack>
#include <vector>
#include "memalign_allocator.h"
#include "vector3.h"
#include "wtime.h"

class nbody_particle{
 private:
  vector3 pos;         // 24-byte
  real mass;           // 32-byte
  vector3 vel;         // 56-byte
  long index;          // 64-byte
  vector3 acc_gravity; // 88-byte
  real phi_gravity;    // 96-byte

 public:
  nbody_particle(){
    assert(sizeof(*this) == 96);
    pos = 0.0;
    vel = 0.0;
    acc_gravity = 0.0;
    phi_gravity = mass =  0.0;
    index = 0;
  }

  static void *operator new[](size_t size){
    return memalign(32, size);
  }
  static void operator delete[](void *p){
    free(p);
  }

  void prefetch(){
    char *addr = (char *)this;
    __builtin_prefetch(addr+0);
    __builtin_prefetch(addr+64);
  }

  void  set_pos(const vector3& new_pos) {pos = new_pos;}
  void  set_vel(const vector3& new_vel) {vel = new_vel;}
  void  set_acc_gravity(const vector3& new_acc) {acc_gravity = new_acc;}
  void  set_phi_gravity(real new_phi) {phi_gravity = new_phi;}
  void  clear_acc_phi_gravity() {acc_gravity = 0.0;phi_gravity = 0.0;}
  void  scale_pos(real scale_factor) {pos *= scale_factor;}
  void  scale_vel(real scale_factor) {vel *= scale_factor;}

  real     get_mass()        const {return mass;}
  vector3  get_pos()         const {return pos;}
  vector3  get_vel()         const {return vel;}
  vector3  get_acc_gravity() const {return acc_gravity;}
  real     get_phi_gravity() const {return phi_gravity;}

  void set_mass(real m) {mass = m;}
  void set_index(int i) {index = i;}
  inline void predict(const real dt){
    const vector3 newvel = vel + (dt*0.5)*acc_gravity;
    pos += dt * newvel;
    vel = newvel;
  }

  void copy(nbody_particle *dst_ptr, bool stream = true) const{
    assert((unsigned long)(dst_ptr) % 16 == 0);
    assert((unsigned long)(this) % 16 == 0);
    typedef float v4sf __attribute__ ((vector_size(16)));
    const v4sf *src = (const v4sf *)this;
    v4sf *dst = (v4sf *)dst_ptr;

    v4sf src0 = src[0];
    v4sf src1 = src[1];
    v4sf src2 = src[2];
    v4sf src3 = src[3];
    v4sf src4 = src[4];
    v4sf src5 = src[5];
    if(!stream){
      dst[0] = src0;
      dst[1] = src1;
      dst[2] = src2;
      dst[3] = src3;
      dst[4] = src4;
      dst[5] = src5;
    }else{
      __builtin_ia32_movntps((float *)(dst + 0), src0);
      __builtin_ia32_movntps((float *)(dst + 1), src1);
      __builtin_ia32_movntps((float *)(dst + 2), src2);
      __builtin_ia32_movntps((float *)(dst + 3), src3);
      __builtin_ia32_movntps((float *)(dst + 4), src4);
      __builtin_ia32_movntps((float *)(dst + 5), src5);
    }
  }

  nbody_particle operator=(const nbody_particle &rhs){
    rhs.copy(this, false);
    return *this;
  }

  void swap(nbody_particle &p){
    typedef float v4sf __attribute__ ((vector_size(16)));
    v4sf *src = (v4sf *)this;
    v4sf *dst = (v4sf *)&p;
    if(src == dst) return;

    v4sf s0 = src[0]; v4sf s1 = src[1];
    v4sf s2 = src[2]; v4sf s3 = src[3];
    v4sf s4 = src[4]; v4sf s5 = src[5];

    v4sf d0 = dst[0]; v4sf d1 = dst[1];
    v4sf d2 = dst[2]; v4sf d3 = dst[3];
    v4sf d4 = dst[4]; v4sf d5 = dst[5];

    src[0] = d0; src[1] = d1;
    src[2] = d2; src[3] = d3;
    src[4] = d4; src[5] = d5;

    dst[0] = s0; dst[1] = s1;
    dst[2] = s2; dst[3] = s3;
    dst[4] = s4; dst[5] = s5;
  }

  void friend accumulate_mutual_gravity(nbody_particle & p1,
                                        nbody_particle & p2,
                                        real eps2);
};

typedef vector3 (nbody_particle::*nbody_VMF_ptr)(void);        // vector3 member function pointer
typedef void (nbody_particle::*nbody_MF_ptr)(const vector3 &); // member function pointer
typedef void (nbody_particle::*nbody_VF_ptr)(void);            // void member function pointer
typedef void (nbody_particle::*nbody_RF_ptr)(real);            // void member function
typedef void (nbody_particle::*nbody_RRF_ptr)(real,real);      // void member function
class bhparticle;
typedef std::vector< nbody_particle, __gnu_cxx::malloc_allocator<nbody_particle> > vec_particle;

class nbody_system
{
 private:
  int nsize;
  nbody_particle *pb;
  double Gflops; 
  vec_particle pbuf1, pbuf2;
 public:
  nbody_particle *get_pbuf1() { return &pbuf1[0]; }
  nbody_particle *get_pbuf2() { return &pbuf2[0]; }
  void swap_pbuf(){
    if(pb == get_pbuf1()){
      pb = get_pbuf2();
    }else if(pb == get_pbuf2()){
      pb = get_pbuf1();
    }else{
      assert(0);
    }
  }
  void set_nsize(const int _nsize);
  int n, _ntot;
  int step;
  real time;
  real timestep;
  real eps2_for_gravity;
  vector3 pos;
  vector3 vel;
  real mass;
  real theta_for_tree;
  int ncrit_for_tree;
  vector3 *xlowp;
  vector3 *xhighp;
  bool do_morton_sort;
  bool do_setup_division;

  nbody_system(){
    static vector3 xlow[NMAXPROC];
    static vector3 xhigh[NMAXPROC];
    xlowp = xlow;
    xhighp = xhigh;
    n = 0;
    step = 0;
    nsize = 0;
    time = 0;
    timestep = 0;
    pb = NULL;
    do_morton_sort = true;
    do_setup_division = true;
  }
 
  void gen_uniform(int);
  void add_particles_to_tree(int nadd);
  void setup_tree();
  void apply_vf(nbody_VF_ptr);
  void apply_vf(nbody_RF_ptr, real);
  void apply_vf(nbody_RRF_ptr, real, real);
  void calculate_gravity();
  void calculate_gravity_direct();
  nbody_particle * get_particle_pointer(){return pb;}
  const nbody_particle * get_particle_pointer() const {return pb;}
  int get_nsize(){return nsize;}
  void morton_sort_using_tree();
  void evaluate_gravity_direct(float eps2, int ntot, bhparticle *bp);

  real rmax(){
    real tmp = 0.0;
    for(int i=0; i<n; i++){
      vector3 pos = pb[i].get_pos();
      tmp = std::max(tmp, pos * pos);
    }
    return sqrt(tmp);
  }
};

#endif
