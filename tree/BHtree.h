/*-----------------------------------------------------------------------------
 *  BHtree : basic class for C++ implementation of BH treecode
 *  J. Makino 1998/12/14
 *  Modified for tree construction through insertion
 *-----------------------------------------------------------------------------
 */

#ifndef BHTREE_H
#define BHTREE_H

#include "nbody_particle.h"

typedef nbody_particle real_particle;
typedef nbody_system real_system;

#define BPWORKMAX 64 

struct float4{
  typedef float  v4sf __attribute__ ((vector_size(16)));
  typedef double v2df __attribute__ ((vector_size(16)));
  static v4sf v4sf_abs(v4sf x){
    typedef int v4si __attribute__ ((vector_size(16)));
    v4si mask = {0x7fffffff, 0x7fffffff, 0x7fffffff, 0x7fffffff};
    return __builtin_ia32_andps(x, (v4sf)mask);
  }
  union{
    v4sf v;
    struct{
      float x, y, z, w;
    };
  };
  float4() : v((v4sf){0.f, 0.f, 0.f, 0.f}) {}
  float4(float x, float y, float z, float w) : v((v4sf){x, y, z, w}) {}
  float4(float x) : v((v4sf){x, x, x, x}) {}
  float4(vector3 vec, real w = 0.0) : v((v4sf){vec[0], vec[1], vec[2], w}) {}
  float4(v4sf _v) : v(_v) {}
  float4 abs(){
    typedef int v4si __attribute__ ((vector_size(16)));
    v4si mask = {0x7fffffff, 0x7fffffff, 0x7fffffff, 0x7fffffff};
    return float4(__builtin_ia32_andps(v, (v4sf)mask));
  }
  operator vector3(){
    return vector3(x, y, z);
  }
  v4sf operator=(const float4 &rhs){
    v = rhs.v;
    return v;
  }
  float4(const float4 &rhs){
    v = rhs.v;
  }
  float get_element(const int dim){
    if(0 == dim) return x;
    if(1 == dim) return y;
    if(2 == dim) return z;
    if(3 == dim) return w;
    // otherwise return NaN
    union{int i; float f;} un = {-1};
    return un.f;
  }
};
#endif

#ifndef PHANTOM
#include "cualloc.h"
typedef std::vector<float4, cuda::cualloc<float4, cudaHostAllocDefault> > vectorfloat4; // pinned memory
#else
typedef std::vector<float4> vectorfloat4; // standared
#endif

extern "C"{
  void force_nwalk(int tid, float4 xi[], float4 xj[], float4 accpot[], int ioff[], int joff[], int nwalk);
  void vforce(int tid, float4 xi[], float4 xj[], float4 accpot[], int ioff[], int joff[], int nwalk);
  void vforce_mp (int pid, int tid, float4 xi[], float4 xj[], float4 accpot[], int ioff[], int joff[], int nwalk);
  void vforce_1st (int pid, int tid, float4 xi[], float4 xj[], int ioff[], int joff[], int nwalk);
  void vforce_2nd (int pid, int tid, float4 accpot[]);
  void vforce_open(int pid, int tid, int num_threads, int num_threads_per_gpu);
  void force_nwalk_sim(int tid, float4 xi[], float4 xj[], float4 accpot[], int ioff[], int joff[], int nwalk);
  void vforce_wait(double total_sec, int nproc, int pid, int tid);
  void vforce_close(int tid);
}

typedef long long BHlong;

class bhparticle
{
private:
  real_particle * rp;
  bhparticle * bhp_next;
  float4 xmcache;
public:
  bhparticle(){
    rp = NULL;
    bhp_next = NULL;
  }
  void set_bhp_next(bhparticle * p){bhp_next = p;}
  bhparticle * get_bhp_next(){return bhp_next ;}
  void set_rp(real_particle * p){
    rp = p;
    bhp_next = NULL;
    xmcache = float4(p->get_pos(), p->get_mass());
//    assert(p->get_mass() > 1.e-10);
//    assert(xmcache.w > 1.e-10);
  }
  void set_rp_addr(real_particle * p){
    rp = p;
  }
  real_particle * get_rp(){return rp ;}
  float4 &get_xmcache(){return xmcache;}
  void set_xmcache(const float4 &v) {xmcache = v; }
  void set_extp(const float4 &v){
    rp = NULL;
    bhp_next = NULL;
    xmcache = v;
//    assert(xmcache.w > 1.e-10);
  }
  vector3 get_pos(){ return vector3(xmcache.x, xmcache.y, xmcache.z); }
};

class bhnode
{
public:
  float4 center;
  float4 cm;
  bhnode     *child;
  bhparticle *bpfirst;
  int nparticle;
  int isleaf;
  int touch;
//  int pad[1];
    
  bhnode(){
    assert(sizeof(*this) == 64);
    center = float4();
    cm = float4();
    child = NULL;
    bpfirst = NULL;
    nparticle = 0;
    isleaf = 1;
    touch = 1;
  }
  void clear(){
    center = float4();
    cm = float4();
    child = NULL;
    bpfirst = NULL;
    nparticle = 0;
    isleaf = 1;
    touch = 1;
  }
  vector3 get_pos(){
    return vector3(center.x, center.y, center.z);
  }
  real get_length(){
    return center.w;
  }
  vector3 get_cmpos(){
    return vector3(cm.x, cm.y, cm.z);
  }
  real get_cmmass(){
    return cm.w;
  }

  void allocate_childrens_and_assign_particles(bhnode *&heap_top, int &heap_remainder, int n_critical);
  void insert_particle(bhparticle *particle, bhnode *&heap_top,int &heap_remainder,int n_critical);
  void assign_root(vector3 root_pos, real length, bhparticle * bp, int nparticle);
  friend int check_and_set_nbl(bhnode * p1,bhnode * p2);
  void set_cm_quantities();
  void set_cm_quantities_parallel();
  void create_essential_tree(vector3 xcen,
			     vector3 size,
			     real theta2,
			     vector3 * &plist,
			     real * &mlist,
			     int & nlist);
  void create_essential_tree_vector(vector3 &center,
				    vector3 &hlen,
				    vector<float4> &xmlist,
				    real theta2);
  int count_real_particles();
  void find_group_list(int ncrit, vector<bhnode*> &group_list);
  void set_i_particles(vector<real_particle *> &iplist,
		       vectorfloat4 &xilist, 
		       float eps2);
#ifndef PHANTOM
  void set_jp(vectorfloat4 &xjlist) const; 
#endif
  void set_jp(std::vector<float4> &xjlist) const; 
  void set_j_particles(vectorfloat4 &xjlist, 
		       bhnode &inode,
		       float theta2);
  void evaluate_gravity_GPU(real theta2, real eps2, int ncrit);
  void morton_dump(real_particle *&pitr);
};

void clear_tree_counters();
bhparticle *get_bp();
void MP_gather_sample_coords(std::vector<float4> &sample_array);
