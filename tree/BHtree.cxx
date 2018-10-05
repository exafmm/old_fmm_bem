/*-----------------------------------------------------------------------------
 *  BHtree : basic class for C++ implementation of BH treecode
 *  J. Makino 1998/12/14
 *  Modified for tree construction through insertion
 *-----------------------------------------------------------------------------
 */

#define NCRIT_FOR_TREE 16

using namespace std;

static int istep;

#define real double
#include "BHtree.h"
#include "nbody.h"
#include "v4sf.h"

static bhparticle *bp = NULL;
bhparticle *get_bp() { return bp; }
static int bhpsize = 0;
static int bnsize = 0;
static bhnode *bn;
static int bnused;
static real total_interactions = 0.0;
static int tree_walks;
static int nisum;

int init_bhp(int nbody,real_particle * rp,int & nbhpsize,bhparticle * &bhp)
{
  int retval = 0;
  if (nbody > nbhpsize || bhp == NULL){
    retval = 1;
    if (bhp != NULL)free(bhp);
    nbhpsize = arraysize_bhp(nbody);
    log_allocate("bhparticle", nbhpsize * sizeof(bhparticle));
    bhp = (bhparticle *)memalign(64, nbhpsize * sizeof(bhparticle));
  }
#pragma omp parallel for
  for(int i = 0; i<nbody; i++){
    bhparticle *p = bhp+i;
    p->set_rp(rp+i);
  }
  return retval;
}

void bhnode::assign_root(vector3 root_pos, real length, bhparticle *bp, int np)
{
  center = float4(root_pos, length);
  bpfirst = bp;
  nparticle = np;
  child = NULL;  
}

void bhnode::set_cm_quantities()
{
  const int nparticle = this->nparticle;
  if(touch == 0) return;
  int i;
  float cmx = 0.0f, cmy = 0.0f, cmz = 0.0f, cmw = 0.0f;
  if (nparticle == 0) return;
  if (isleaf){
    bhparticle *bp = bpfirst;
    for(i = 0; i < nparticle; i++){
      const float4 &xm = bp->get_xmcache();
      cmx += xm.w * xm.x;
      cmy += xm.w * xm.y;
      cmz += xm.w * xm.z;
      cmw += xm.w;
      bp=bp->get_bhp_next(); 
    }
  }else{
    if (child != NULL){
      for(i=0;i<8;i++){
        child[i].set_cm_quantities();
      }
      for(i=0;i<8;i++){
        const float4 &xm = child[i].cm;
        cmx += xm.w * xm.x;
        cmy += xm.w * xm.y;
        cmz += xm.w * xm.z;
        cmw += xm.w;
      }
    }
  }
  cmx *= 1.f/cmw;
  cmy *= 1.f/cmw;
  cmz *= 1.f/cmw;
  cm.x = cmx;
  cm.y = cmy;
  cm.z = cmz;
  cm.w = cmw;
  touch = 0;
}

bhnode *get_bhroot()
{
  return bn;
}

void set_cm_quantities_for_default_tree()
{
  bn->set_cm_quantities();
}

static inline int octant_index(const float4 &center, const float4 &pos){
  int i = __builtin_ia32_movmskps(
    (float4::v4sf)__builtin_ia32_cmpltps(center.v, pos.v));
  return 7 & i;
}

void bhnode::allocate_childrens_and_assign_particles(
  bhnode *&heap_top, 
  int &heap_remainder, 
  int n_critical)
{
  const int nparticle = this->nparticle;
  static bhparticle *bhparray[BPWORKMAX];
  assert(nparticle <= BPWORKMAX);
  bhparticle *p = bpfirst;
  for(int i=0; i<nparticle; i++){
    bhparray[i] = p;
    p = p->get_bhp_next();
  }
  int nbhp = nparticle;
  assert(child == NULL);
  child = heap_top;
  heap_top +=8;
  heap_remainder -=8 ;
  assert(heap_remainder >= 0);

  const float l = center.w;
  const float4::v4sf lvec = {l, l, l, l};
  static const float4::v4sf table[8] = {
    {-0.25, -0.25, -0.25, -0.5},
    {+0.25, -0.25, -0.25, -0.5},
    {-0.25, +0.25, -0.25, -0.5},
    {+0.25, +0.25, -0.25, -0.5},
    {-0.25, -0.25, +0.25, -0.5},
    {+0.25, -0.25, +0.25, -0.5},
    {-0.25, +0.25, +0.25, -0.5},
    {+0.25, +0.25, +0.25, -0.5},
  };
  for(int ic=0; ic<8; ic++){
    bhnode node;
    node.center.v = center.v + table[ic] * lvec;
    child[ic] = node;
  }
  for(int i=0; i<nbhp; i++){
    bhparticle *p = bhparray[i];
    int ic = octant_index(center, p->get_xmcache());
    p->set_bhp_next(child[ic].bpfirst);
    child[ic].bpfirst = p;
    child[ic].nparticle++;
    child[ic].touch = 1;
  }
  for(int ic=0; ic<8; ic++){
    if(child[ic].nparticle >= n_critical){
      child[ic].allocate_childrens_and_assign_particles(heap_top, heap_remainder, n_critical);
    }
  }
  isleaf = 0;
}

void bhnode::insert_particle(bhparticle * particle,
           bhnode * & heap_top,
           int & heap_remainder,
           int n_critical)
{
  bhnode * current = this;
  float4 xm = particle->get_xmcache();
  for(;;){
    if(current->child){ // is not leaf
      int ic = octant_index(current->center, xm);
      bhnode *bhnext = &current->child[ic];
      current->nparticle++;
      current->touch = 1;
      current = bhnext;
      continue;
    }
    assert(current->isleaf);
    if(current->nparticle+1 < n_critical){
      particle->set_bhp_next(current->bpfirst);
      current->bpfirst = particle;
      current->nparticle++;
      current->touch = 1;
      return;
    }else{
      particle->set_bhp_next(current->bpfirst);
      current->bpfirst = particle;
      current->nparticle++;
      current->touch = 1;
      current->allocate_childrens_and_assign_particles(heap_top,heap_remainder, n_critical);
      return;
    }
  }
}
    
void real_system::add_particles_to_tree(int nadd)
{
  bhnode *btmp = bn+bnused;
  int heap_remainder = bnsize-bnused;
  const int nbhpsize = arraysize_bhp(n);
  assert(nbhpsize >= n+nadd);
  for(int i=n; i<n+nadd; i++){
    bn->insert_particle(bp+i, btmp, heap_remainder, NCRIT_FOR_TREE);
  }
  bnused = bnsize - heap_remainder; 
}

void real_system::setup_tree()
{
  ::istep = this->step;
  real rsize = 1.0;
  init_bhp(n, get_particle_pointer(), bhpsize, bp);
  real rmax = -xlowp[0][0];
  while(rsize < rmax) rsize*= 2;
    
  int expected_bnsize =  (int)(bhpsize*0.25 + 100);
  if (bnsize < int(0.9 * expected_bnsize)){
    if (bnsize != 0) free(bn);
    bnsize = expected_bnsize;
    log_allocate("bhnode", bnsize * sizeof(bhnode));
    bn = (bhnode *)memalign(64, bnsize * sizeof(bhnode));
  }
  for(int j=0; j<bnsize; j++) bn[j].clear();
  bn->assign_root(vector3(0.0), rsize*2, NULL, 0);
  bhnode * btmp = bn + 1;
  int heap_remainder = bnsize-1;
  for(int i = 0;i<n;i++){
    bn->insert_particle(bp+i, btmp, heap_remainder, NCRIT_FOR_TREE);
  }
  bnused = bnsize - heap_remainder; 
}

void real_system::morton_sort_using_tree(){
  nbody_particle *pitr = pb;
  bn->morton_dump(pitr);
}

extern double get_ninteractions() {
  return total_interactions;
}

void clear_tree_counters()
{
  total_interactions = 0;
  tree_walks = 0;
  nisum = 0;
}

static void calc_opening_status_let(
  int opening_status[8], 
  const float4 &center, 
  const float4 &hlen, 
  const bhnode child[8], 
  const float4 &theta2){
  v4sf_x4 ri(center.v);
  v4sf_x4 hi(hlen.v);
  v4sf_x4 rj0(child[0].cm.v, child[1].cm.v, child[2].cm.v, child[3].cm.v, true);
  rj0.w = v4sf(child[0].center.w).val;
  v4sf_x4 rj1(child[4].cm.v, child[5].cm.v, child[6].cm.v, child[7].cm.v, true);
  rj1.w = v4sf(child[0].center.w).val;

  v4sf os0, os1;
  {
    const v4sf_x4 &rj = rj0;
    v4sf dx = v4sf(rj.x - ri.x).abs();
    v4sf dy = v4sf(rj.y - ri.y).abs();
    v4sf dz = v4sf(rj.z - ri.z).abs();
    v4sf l2 = rj.w * rj.w;
    dx -= hi.x;
    dy -= hi.y;
    dz -= hi.z;
    v4sf dx2 = dx + dx.abs();
    v4sf dy2 = dy + dy.abs();
    v4sf dz2 = dz + dz.abs();
    v4sf r2 = dx2*dx2 + dy2*dy2 + dz2*dz2;
    v4sf tmp = (r2 * (v4sf(0.25f) * theta2.v)) > l2;
    v4sf hj = v4sf(0.5f) * rj.w;
    v4sf status = tmp & ( (dx >= hj) | (dy >= hj) | (dz >= hj) );
    os0 = status;
  }
  {
    const v4sf_x4 &rj = rj1;
    v4sf dx = v4sf(rj.x - ri.x).abs();
    v4sf dy = v4sf(rj.y - ri.y).abs();
    v4sf dz = v4sf(rj.z - ri.z).abs();
    v4sf l2 = rj.w * rj.w;
    dx -= hi.x;
    dy -= hi.y;
    dz -= hi.z;
    v4sf dx2 = dx + dx.abs();
    v4sf dy2 = dy + dy.abs();
    v4sf dz2 = dz + dz.abs();
    v4sf r2 = dx2*dx2 + dy2*dy2 + dz2*dz2;
    v4sf tmp = (r2 * (v4sf(0.25f) * theta2.v)) > l2;
    v4sf hj = v4sf(0.5f) * rj.w;
    v4sf status = tmp & ( (dx >= hj) | (dy >= hj) | (dz >= hj) );
    os1 = status;
  }
  *(v4sf *)(opening_status+0) = os0;
  *(v4sf *)(opening_status+4) = os1;
}

void bhnode::set_jp(std::vector<float4> &xjlist) const
{       
  if(isleaf){
    bhparticle * bp = bpfirst;
    for(int i=0; i<nparticle; i++){
      xjlist.push_back(bp->get_xmcache());
      bp = bp->get_bhp_next();
    }
  }else{  
    if(child != NULL){
      for(int ic=0; ic<8; ic++){
        child[ic].set_jp(xjlist);
      }
    }
  }
}

void bhnode::create_essential_tree_vector(
  vector3 &_center,
  vector3 &_hlen,
  vector<float4> &xmlist,
  real _theta2)
{
  if(nparticle == 0) return;
  float4 center(_center);
  float4 hlen  (_hlen);
  float4 theta2 = (float(_theta2));

  std::stack <const bhnode *> stack;
  assert(!this->isleaf);
  stack.push(this);
  while(!stack.empty()){
    const bhnode *node = stack.top();
    stack.pop();
    int opening_status[8] __attribute__ ((aligned(16)));
    calc_opening_status_let(opening_status, center, hlen, node->child, theta2);
    for(int ic=0; ic<8; ic++){
      const bhnode &cnode = node->child[ic];
      if(cnode.nparticle == 0) continue;
      if(opening_status[ic]){
        xmlist.push_back(cnode.cm);
      }else{
        if(cnode.isleaf){
          cnode.set_jp(xmlist);
        }else{
          stack.push(&cnode);
        }
      }
    }
  }
}

void evaluate_gravity_using_default_tree_and_list(real theta2, real eps2, int ncrit)
{
  try{
    bn->evaluate_gravity_GPU(theta2, eps2, ncrit);
  }catch(bad_alloc){
    int myid   = MP_myprocid();
    cout << myid << " : bad alloc 2" << endl;
    MP_Abort(-1);
  }
}

int bhnode::count_real_particles(){
  int ni=0;
  if(isleaf){
    bhparticle * bp = bpfirst;
    while(bp != NULL){
      real_particle * p =  bp->get_rp();
      if(p != NULL){ // An external particle has NULL
        ni++;
      }
      bp = bp->get_bhp_next();
    }
    return ni;
  }else{
    if(child != NULL){
      for(int ic=0;ic<8;ic++){
        ni += child[ic].count_real_particles();
      }
    }
  }
  return ni;
}

void bhnode::find_group_list(int ncrit, vector<bhnode*> &group_list){
  if((nparticle > ncrit) && (!isleaf)){
    if(child != NULL){
      for(int ic=0;ic<8;ic++){
        child[ic].find_group_list(ncrit, group_list);
      }
    }
  }else{
    if(count_real_particles() > 0){
      group_list.push_back(this);
    }
  }
}

void bhnode::set_i_particles(
    vector<real_particle *> &iplist,
    vectorfloat4 &xilist, 
    float eps2){
  if(isleaf){
    bhparticle * bp = bpfirst;
    while(bp != NULL){
      real_particle * p =  bp->get_rp();
      if(p != NULL){ // An external particle has NULL
        iplist.push_back(p);
        float4 ip = bp->get_xmcache();
        ip.w = eps2;
        xilist.push_back(ip);
      }
      bp = bp->get_bhp_next();
    }
  }else{
    if(child != NULL){
      for(int ic=0;ic<8;ic++){
        child[ic].set_i_particles(iplist, xilist, eps2);
      }
    }
  }
}

static void calc_opening_status(int opening_status[8], const float4 &ipos, const bhnode child[8], const float _theta2){
  v4sf theta2 = _theta2;
  v4sf_x4 ri(ipos.v);
  v4sf_x4 rj0(child[0].cm.v, child[1].cm.v, child[2].cm.v, child[3].cm.v, true);
  rj0.w = v4sf(child[0].center.w).val;
  v4sf_x4 rj1(child[4].cm.v, child[5].cm.v, child[6].cm.v, child[7].cm.v, true);
  rj1.w = v4sf(child[0].center.w).val;

  v4sf os0, os1;
  {
    const v4sf_x4 &rj = rj0;
    v4sf dx = v4sf(rj.x - ri.x).abs();
    v4sf dy = v4sf(rj.y - ri.y).abs();
    v4sf dz = v4sf(rj.z - ri.z).abs();
    v4sf hi = v4sf(0.5f) * ri.w;
    v4sf hj = v4sf(0.5f) * rj.w;
    v4sf l2 = rj.w * rj.w;
    dx -= hi;
    dy -= hi;
    dz -= hi;
    v4sf dx2 = dx + dx.abs();
    v4sf dy2 = dy + dy.abs();
    v4sf dz2 = dz + dz.abs();
    v4sf r2 = dx2*dx2 + dy2*dy2 + dz2*dz2;
    v4sf tmp = (r2 * (v4sf(0.25f) * theta2)) > l2;
    v4sf status = tmp & ( (dx >= hj) | (dy >= hj) | (dz >= hj) );
    os0 = status;
  }
  {
    const v4sf_x4 &rj = rj1;
    v4sf dx = v4sf(rj.x - ri.x).abs();
    v4sf dy = v4sf(rj.y - ri.y).abs();
    v4sf dz = v4sf(rj.z - ri.z).abs();
    v4sf hi = v4sf(0.5f) * ri.w;
    v4sf hj = v4sf(0.5f) * rj.w;
    v4sf l2 = rj.w * rj.w;
    dx -= hi;
    dy -= hi;
    dz -= hi;
    v4sf dx2 = dx + dx.abs();
    v4sf dy2 = dy + dy.abs();
    v4sf dz2 = dz + dz.abs();
    v4sf r2 = dx2*dx2 + dy2*dy2 + dz2*dz2;
    v4sf tmp = (r2 * (v4sf(0.25f) * theta2)) > l2;
    v4sf status = tmp & ( (dx >= hj) | (dy >= hj) | (dz >= hj) );
    os1 = status;
  }
  *(v4sf *)(opening_status+0) = os0;
  *(v4sf *)(opening_status+4) = os1;
}

#ifndef PHANTOM
void bhnode::set_jp(vectorfloat4 &xjlist) const{
  if(isleaf){
    bhparticle * bp = bpfirst;
    for(int i=0; i<nparticle; i++){
      xjlist.push_back(bp->get_xmcache());
      bp = bp->get_bhp_next();
    }
  }else{
    if(child != NULL){
      for(int ic=0;ic<8;ic++){
        child[ic].set_jp(xjlist);
      }
    }
  }
}  
#endif

void bhnode::set_j_particles(vectorfloat4 &xjlist, bhnode &inode, float theta2){
  std::stack <const bhnode *> stack;
  stack.push(this);
  while(!stack.empty()){
    const bhnode *node = stack.top();
    stack.pop();
    int opening_status[8] __attribute__ ((aligned(16)));
    calc_opening_status(opening_status, inode.center, node->child, theta2);
    for(int ic=0; ic<8; ic++){
      const bhnode &cnode = node->child[ic];
      if(cnode.nparticle == 0) continue;
      if(opening_status[ic]){
        xjlist.push_back(cnode.cm);
      }else{
        if(cnode.isleaf || (&cnode == &inode)){
          cnode.set_jp(xjlist);
        }else{
          stack.push(&cnode);
        }
      }
    }
  }
}

static void assign_gpu_force(
    vector<real_particle *> &iplist,
    const vectorfloat4 &accplist,
    float4 &amax)
{
  for(int i=0; i<int(iplist.size()); i++){
    real_particle *p = iplist[i];
    if(p){
      const float4 &accp = accplist[i];
      typedef double v2df __attribute__ ((vector_size(16)));
      v2df a0 = {accp.x, accp.y};
      v2df a1 = {accp.z, accp.w};
      __builtin_ia32_movntpd((double *)p +  8, a0);
      __builtin_ia32_movntpd((double *)p + 10, a1);
      amax.v = __builtin_ia32_maxps(amax.v,  accp.v*accp.v);
    }
  }
}

const int MAX_WALK = 256;
const int NIALLOC =  2000 * MAX_WALK;
const int NJALLOC = 20000 * MAX_WALK;
#include "omp_lock.h"
#ifdef PHANTOM // sync version
void bhnode::evaluate_gravity_GPU(real theta2, real eps2, int ncrit){
  const int TMAX = 4;
  static vector<bhnode*> group_list;
  static vector<real_particle *> iplist[TMAX];
  static vectorfloat4 xilist[TMAX]; 
  static vectorfloat4 xjlist[TMAX]; 
  static vectorfloat4 accplist[TMAX]; 
  static bool firstcall = true;
  float4 amax(0.f);
  if(firstcall){
#pragma omp parallel
    {
      int total_proc_count = MP_proccount();
      int pid = MP_myprocid();
      int tid = omp_get_thread_num();

      vforce_wait(total_proc_count*3.5/32.0, total_proc_count, pid, tid);
      vforce_open(pid, tid, omp_get_num_threads(), 2); // for normal cudaMalloc
      iplist[tid].reserve(1<<16);
      xilist[tid].reserve(1<<16);
      xjlist[tid].reserve(1<<16);
      accplist[tid].reserve(1<<16);
    }
    firstcall = false;
  }
  const int IALIGN   = 1;
  const int JALIGN   = 4;
  static int ioff[TMAX][MAX_WALK+1];
  static int joff[TMAX][MAX_WALK+1];
  group_list.clear();
  find_group_list(ncrit, group_list);
  int ng = group_list.size();
  double thost = 0.0, tgpu = 0.0;
  omp_lock lock;

#pragma omp parallel for schedule(dynamic) reduction(+ : thost, tgpu, total_interactions, tree_walks, nisum) 
  for(int ig=0; ig<ng; ig+=MAX_WALK){
    int tid = omp_get_thread_num();
    double t0 = get_wtime();
    int nwalk = min(MAX_WALK, ng-ig);
    iplist[tid].clear();
    xilist[tid].clear();
    xjlist[tid].clear();
    ioff[tid][0] = 0;
    joff[tid][0] = 0;
    for(int iw=0; iw<nwalk; iw++){
      bhnode *inode = group_list[ig + iw];
      inode->set_i_particles(iplist[tid], xilist[tid], eps2);
      while(xilist[tid].size() % IALIGN){
        iplist[tid].push_back(NULL);
        xilist[tid].push_back(float4(0.f, 0.f, 0.f, 0.f));
      }
      bn->set_j_particles(xjlist[tid], *inode, theta2);
      int ni = xilist[tid].size() - ioff[tid][iw];
      int nj = xjlist[tid].size() - joff[tid][iw];
      while(xjlist[tid].size() % JALIGN){ 
        xjlist[tid].push_back(float4(0.f, 0.f, 0.f, 0.f));
      }
      ioff[tid][iw+1] = xilist[tid].size();
      joff[tid][iw+1] = xjlist[tid].size();
      total_interactions += ni * nj;
      tree_walks++;
      nisum += ni;
    }
    xilist  [tid].resize(xilist[tid].size() +4);//xxx
    accplist[tid].resize(xilist[tid].size());
    double t1 = get_wtime();
    lock.set_lock();
    vforce_mp(MP_myprocid(), tid, &xilist[tid][0], &xjlist[tid][0], &accplist[tid][0], ioff[tid], joff[tid], nwalk);

    lock.unset_lock();
    double t2 = get_wtime();

    assign_gpu_force(iplist[tid], accplist[tid], amax);

    thost += t1 - t0;
    tgpu  += t2 - t1;
  } // end omp paralles for
}
#  else // async version
void bhnode::evaluate_gravity_GPU(real theta2, real eps2, int ncrit){
  const int TMAX = 8;
  static vector<bhnode*> group_list;
  static vector<real_particle *> iplist[TMAX];
  static vector<real_particle *> iplist_prev[TMAX];
  static vectorfloat4 xilist[TMAX]; 
  static vectorfloat4 xjlist[TMAX]; 
  static vectorfloat4 accplist[TMAX]; 
  static bool firstcall = true;
  static int nthreads;
  float4 amax(0.f);
  omp_lock lock;
  if(firstcall){
#pragma omp parallel
    {
      int total_proc_count = MP_proccount();
      int pid = MP_myprocid();
      int tid = omp_get_thread_num();

      vforce_wait(total_proc_count*3.5/512.0, total_proc_count, pid, tid); // wait for 0.11 second for each open

      lock.set_lock();
      {
        const int nthread_per_gpu = 1;
        vforce_open(pid, tid, omp_get_num_threads(), nthread_per_gpu); // for normal cudaMalloc
        const size_t isize = NIALLOC;
        const size_t jsize = NJALLOC;
        iplist     [tid].reserve(isize);
        iplist_prev[tid].reserve(isize);
        xilist     [tid].reserve(isize);
        xjlist     [tid].reserve(jsize);
        accplist   [tid].reserve(isize);
      }
      lock.unset_lock();
#pragma omp master
      nthreads = omp_get_num_threads();
    }
    firstcall = false;
    MP_sync();
  }
  const int IALIGN   = 1;
  const int JALIGN   = 128;
  static int ioff[TMAX][MAX_WALK+1];
  static int joff[TMAX][MAX_WALK+1];
  group_list.clear();
  find_group_list(ncrit, group_list);
  int ng = group_list.size();
  double thost = 0.0, tgpu = 0.0;

  static int goff[TMAX+1];
  for(int i=0; i<=nthreads; i++){
    goff[i] = (i * ng) / nthreads;
  }
#pragma omp parallel reduction(+ : thost, tgpu, total_interactions, tree_walks, nisum) 
  {
    int tid = omp_get_thread_num();
    int gstart = goff[tid];
    int gend   = goff[tid+1];

    for(int ig=gstart; ig<gend; ig+= MAX_WALK){
    int nwalk = min(MAX_WALK, gend-ig);

    double t0 = get_wtime();
    iplist[tid].clear();
    xilist[tid].clear();
    xjlist[tid].clear();
    ioff[tid][0] = 0;
    joff[tid][0] = 0;
    for(int iw=0; iw<nwalk; iw++){
      bhnode *inode = group_list[ig + iw];
      inode->set_i_particles(iplist[tid], xilist[tid], eps2);
      while(xilist[tid].size() % IALIGN){
        iplist[tid].push_back(NULL);
        xilist[tid].push_back(float4(0.f, 0.f, 0.f, 0.f));
      }
      bn->set_j_particles(xjlist[tid], *inode, theta2);
      int ni = xilist[tid].size() - ioff[tid][iw];
      int nj = xjlist[tid].size() - joff[tid][iw];
      while(xjlist[tid].size() % JALIGN){ 
        xjlist[tid].push_back(float4(0.f, 0.f, 0.f, 0.f));
      }
      ioff[tid][iw+1] = xilist[tid].size();
      joff[tid][iw+1] = xjlist[tid].size();
      total_interactions += ni * nj;
      tree_walks++;
      nisum += ni;
    }
    xilist [tid].resize(xilist[tid].size() +4);

    double t1 = get_wtime();

    if(!iplist_prev[tid].empty()){
      vforce_2nd(MP_myprocid(), tid, &accplist[tid][0]);
      assign_gpu_force(iplist_prev[tid], accplist[tid],  amax);
      iplist_prev[tid].clear();
    }
    vforce_1st(MP_myprocid(), tid, &xilist[tid][0], &xjlist[tid][0], ioff[tid], joff[tid], nwalk);
    accplist[tid].resize(xilist[tid].size());
    iplist[tid].swap(iplist_prev[tid]);

    double t2 = get_wtime();

    thost += t1 - t0;
    tgpu  += t2 - t1;
    }
  }
#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    if(!iplist_prev[tid].empty()){
      vforce_2nd(MP_myprocid(), tid, &accplist[tid][0]);
      assign_gpu_force(iplist_prev[tid], accplist[tid],  amax);
      iplist_prev[tid].clear();
    }
  }
  if(__builtin_isnanf(amax.x + amax.y + amax.z + amax.w)){
    static int counter = 0;
    int myid = MP_myprocid();
    const char *myname = MP_get_hostname();
    printf("!!! %d WARING !!!!: NaN is Detected, retrying evaluate_gravity_GPU (hostname = %s, counter = %d)\n", myid, myname, ++counter);
    clear_tree_counters();
    evaluate_gravity_GPU(theta2, eps2, ncrit);
  }
  if(sqrtf(amax.x + amax.y + amax.z) > 1000.){
    static int counter = 0;
    int myid = MP_myprocid();
    const char *myname = MP_get_hostname();
    printf("!!! %d WARING !!!!: Irregular force is Detected, retrying evaluate_gravity_GPU (hostname = %s, counter = %d)\n", myid, myname, ++counter);
    clear_tree_counters();
    evaluate_gravity_GPU(theta2, eps2, ncrit);
  }
}
#  endif

void bhnode::morton_dump(real_particle *&pitr){
  const int nparticle = this->nparticle;
  if(isleaf){
    static bhparticle *bhparray[BPWORKMAX];
    bhparticle *p = bpfirst;
    for(int i=0; i<nparticle; i++){
      bhparray[i] = p;
      p->get_rp()->prefetch();
      p = p->get_bhp_next();
    }
    for(int i=0; i<nparticle; i++){
      bhparticle *bp = bhparray[i];
      real_particle *p =  bp->get_rp();
      p->copy(pitr); // store the particle
      bp->set_rp_addr(pitr);
      pitr++;
    }
  }else{
    if(child != NULL){
      for(int ic=0;ic<8;ic++){
        child[ic].morton_dump(pitr);
      }
    }
  }
}
