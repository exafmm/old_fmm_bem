// Copyright(C) 2010 by 
// Tsuyoshi Hamada <hamada@progrape.jp>
// Keigo Nitadori <nitadori@margaux.astron.s.u-tokyo.ac.jp>
// Rio Yokota <rio.yokota@bristol.ac.uk>

using namespace std;

#define real double
#include "BHtree.h"
#include "nbody.h"
#include <omp.h>

typedef nbody_particle real_particle;
typedef nbody_system real_system;
void set_cm_quantities_for_default_tree();
void clear_tree_counters();

void evaluate_gravity_using_default_tree_and_list(real theta2,
						  real eps2,
						  int ncrit);

void real_system::evaluate_gravity_direct(
					  float eps2, 
					  int ntot, 
					  bhparticle *bp){
  int ioff[] = {0, n};
  int joff[] = {0, ntot};
  vectorfloat4 xilist;
  vectorfloat4 xjlist;
  vectorfloat4 accplist;
  xilist.reserve(n);
  xjlist.reserve(ntot);
  accplist.resize(n);
  for(int i=0; i<n; i++){
    xilist.push_back(float4(pb[i].get_pos(), eps2));
  }
  for(int j=0; j<ntot; j++){
    xjlist.push_back(bp[j].get_xmcache());
  }
  vforce_mp(MP_myprocid(), 0, &xilist[0], &xjlist[0], &accplist[0], ioff, joff, 1);
  for(int i=0; i<n; i++){
    float4 &accp = accplist[i];
    pb[i].set_acc_gravity(vector3(accp.x, accp.y, accp.z));
    pb[i].set_phi_gravity(accp.w);
  }
}

static int num_threads(){
#ifdef _OPENMP
  int nth;
#pragma omp parallel
  {
#pragma omp master
    nth = omp_get_num_threads();
  }
  return nth;
#else
  return 1;
#endif
}

void real_system::calculate_gravity()
{
  const int nthread = num_threads();
  try{
    setup_tree();
  }catch(bad_alloc){
    int myid   = MP_myprocid();
    cout << myid << " : bad alloc 3" << endl;
    MP_Abort(-1);
  }
  if(nthread == 1){
    if(do_morton_sort){
      swap_pbuf();
      morton_sort_using_tree();
    }
  }
  MP_sync();
  set_cm_quantities_for_default_tree();
  MP_sync();
  int ntot;
  real theta2 = theta_for_tree * theta_for_tree;
  if(nthread > 1){
#pragma omp parallel
    {
      const int tid = omp_get_thread_num();
      if(0 == tid){
        exchange_local_essential_trees_alltoall_detune(
            pb, get_bp(), n, nsize,theta2, xlowp, xhighp, get_bhroot(), ntot);
      }
      if(1 == tid){
        swap_pbuf();
        morton_sort_using_tree();
      }
    }
  }else{
    int nbhpmax = arraysize_bhp(n);
    exchange_local_essential_trees_alltoall_detune(
        pb, get_bp(), n, nbhpmax,theta2, xlowp, xhighp, get_bhroot(), ntot);
  }
  _ntot = ntot;
  if(ntot > n){
    add_particles_to_tree(ntot-n);
    set_cm_quantities_for_default_tree();
  }
  MP_sync();
  clear_tree_counters();
  evaluate_gravity_using_default_tree_and_list(theta_for_tree*theta_for_tree, eps2_for_gravity,ncrit_for_tree);
#if 0
  double err,ern=0,erd=0,pbd[10000];
  for(int i=0; i<this->n; i++) pbd[i] = pb[i].get_acc_gravity()[0];
  bhparticle *get_bp();
  evaluate_gravity_direct(eps2_for_gravity,ncrit_for_tree,get_bp());
  for(int i=0; i<this->n; i++) {
    erd += (pb[i].get_acc_gravity()[0]-pbd[i])*(pb[i].get_acc_gravity()[0]-pbd[i]);
    ern += pbd[i]*pbd[i];
  }
  err = sqrt(erd/ern);
  cout << "L2 norm error : " << err << endl;
#endif
}

void real_system::calculate_gravity_direct()
{
  apply_vf(&real_particle::clear_acc_phi_gravity);
  int i, j;
  real_particle * pi;
  real_particle * pj;
  for(i = 0,  pi = &(pb[0]); i<n-1; i++,pi++){
    for(j = i+1,  pj = pi+1; j<n; j++,pj++){
      accumulate_mutual_gravity(*pi, *pj, eps2_for_gravity);
    }
  }
}

void accumulate_mutual_gravity(real_particle & p1, real_particle & p2, real eps2)
{
  vector3 dx = p1.pos-p2.pos;
  double r2inv = 1/(dx*dx+eps2);
  double rinv  = sqrt(r2inv);
  double r3inv = r2inv*rinv;
  p1.phi_gravity -= p2.mass*rinv;
  p2.phi_gravity -= p1.mass*rinv;
  p1.acc_gravity -= p2.mass*r3inv*dx;
  p2.acc_gravity += p1.mass*r3inv*dx;
}
