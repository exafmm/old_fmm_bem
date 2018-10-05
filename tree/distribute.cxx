//
// distribute.C
//
// Version 2001/9/4 J. Makino
//

using namespace std;
#define real double
#include "BHtree.h"
#include "nbody.h"
#include <boost/format.hpp>

typedef nbody_particle real_particle;
typedef nbody_system real_system;

static void create_division(const int n, int &nx, int &ny, int &nz)
{ 
  int n0, n1; 
  n0 = (int)pow(n+0.1,0.33333333333333333333);
  while(n%n0)n0--;
  nx = n0;
  n1 = n/nx;
  n0 = (int)sqrt(n1+0.1);
  while(n1%n0)n0++;
  ny = n0; nz = n1/n0;
  if (nz > ny){
    std::swap(ny, nz);
  }
  if (ny > nx){
    std::swap(nx, ny);
  }
  if (nz > ny){
    std::swap(ny, nz);
  }
// IB topology dependent tuning
  if(96 == n){
    nx = 6;
    ny = 4;
    nz = 4;
  }
  if(192 == n){
    nx = 8;
    ny = 6;
    nz = 4;
  }
  if(288 == n){
    nx = 6;
    ny = 6;
    nz = 8;
  }
  if(384 == n){
    nx = 8;
    ny = 6;
    nz = 8;
  }
  if(480 == n){
    nx = 10;
    ny = 6;
    nz = 8;
  }
  if(576 == n){
    nx = 12;
    ny = 6;
    nz = 8;
  }
  if (nx*ny*nz != n){
    cerr << "create_division: Intenal Error " << n << " " << nx
         << " " << ny << " " << nz <<endl;
    MP_Abort(1);
  }
}

real GLOBAL_ramx = 0.0;

template <int mask> struct CmpFloat4{
  bool operator()(const float4 &lhs, const float4 &rhs){
    return 
    mask & __builtin_ia32_movmskps(
      (float4::v4sf)__builtin_ia32_cmpltps(lhs.v, rhs.v));
  }
};

struct RecursiveDivision{
  const int dim;  // 0, 1 or 2
  const int ndiv;
  std::vector<float> boundary;
  std::vector<RecursiveDivision> child;
  std::vector<float4 *> new_begin;

  RecursiveDivision(const int _dim, const int * const ndivp) : dim(_dim), ndiv(*ndivp)
  {
    boundary.resize(ndiv + 1);
    if(dim < 2){
      child.resize(ndiv, RecursiveDivision(dim+1, ndivp+1));
      new_begin.resize(ndiv + 1);
    }
  }

  ~RecursiveDivision(){
  }

  void operator=(const RecursiveDivision &rd){
    *(const_cast<int *>(&dim) ) = rd.dim;
    *(const_cast<int *>(&ndiv)) = rd.ndiv;
  }

  void calculate_boundary(float4 * const begin, float4 * const end, const float rmin[], const float rmax[])
  {
    if(0==dim) std::sort(begin, end, CmpFloat4<1>());
    if(1==dim) std::sort(begin, end, CmpFloat4<2>());
    if(2==dim) std::sort(begin, end, CmpFloat4<4>());
    
    float wtot = 0.f;
    for(float4 *it = begin; it != end; ++it){
      wtot += (*it).w;
    }
    float wsum = 0.0f;
    float4 *it = begin;
    for(int ix=0; ix<ndiv; ix++){
      const float wc = wtot * float(ix)/float(ndiv);
      while(wsum < wc){
        wsum += (*it++).w;
      }
      if(ix > 0){
        assert(it-2 >= begin); // at least 2 samples are needed to evaluate boundary[1];
        const float wb = wsum;
        const float wa = wsum - it[-1].w;
        assert((wa <= wc) && (wc <= wb));
        const float xa = it[-2].get_element(dim);
        const float xb = it[-1].get_element(dim);
        const float xc = ((wb - wc) * xa + (wc - wa) * xb) / (wb - wa);

        boundary[ix] = xc;
        if(dim < 2){
          new_begin[ix] = it;
        }
      }
    }
    boundary[0]    = rmin[dim];
    boundary[ndiv] = rmax[dim];
    if(dim < 2){
      new_begin[0]    = begin;
      new_begin[ndiv] = end;
      // recursive call
      for(int ix=0; ix<ndiv; ix++){
        child[ix].calculate_boundary(new_begin[ix], new_begin[ix+1], rmin, rmax);
      }
    }
  }
  void write_boundary(const int ibox[], double xlow[], double xhigh[])
  {
    *xlow  = boundary[0 + *ibox];
    *xhigh = boundary[1 + *ibox];
    if(dim < 2){
      child[*ibox].write_boundary(ibox+1, xlow+1, xhigh+1);
    }
  }
};

static void determine_division_with_weight( // nitadori's version
    const int np, // number of particles
    std::vector<float4> &pos, // positions of particles
    const int nx,
    const int ny,
    const int nz,
    const real rmax,
    vector3 xlow[], // left-bottom coordinate of divisions
    vector3 xhigh[]) // size of divisions
{

  const int ndiv[3] = {nx, ny, nz};
  RecursiveDivision rd(0, ndiv);
  const float xmin[3] = {-rmax, -rmax, -rmax};
  const float xmax[3] = {+rmax, +rmax, +rmax};
  rd.calculate_boundary(&pos[0], &pos[np], xmin, xmax);
  for(int ii=0, ix=0; ix<nx; ix++){
    for(int iy=0; iy<ny; iy++){
      for(int iz=0; iz<nz; iz++, ii++){
        const int ibox[3] = {ix, iy, iz};
        rd.write_boundary(ibox, &xlow[ii][0], &xhigh[ii][0]);
      }
    }
  }
}

inline bool isinbox(const vector3 &pos, const vector3 &xlow, const vector3 &xhigh)
{
  return !(
  (pos[0] < xlow [0]) || 
  (pos[1] < xlow [1]) || 
  (pos[2] < xlow [2]) || 
  (pos[0] > xhigh[0]) || 
  (pos[1] > xhigh[1]) || 
  (pos[2] > xhigh[2]) );
}
struct Boundary{
  double xlow, xhigh;
  double ylow, yhigh;
  double zlow, zhigh;
  Boundary(const vector3 &low, const vector3 &high){
    xlow = low[0]; xhigh = high[0];
    ylow = low[1]; yhigh = high[1];
    zlow = low[2]; zhigh = high[2];
  }
  bool isinbox(const vector3 &pos) const {
    return !(
      (pos[0] < xlow ) || 
      (pos[1] < ylow ) || 
      (pos[2] < zlow ) || 
      (pos[0] > xhigh) || 
      (pos[1] > yhigh) || 
      (pos[2] > zhigh) );
  }
};

static void collect_sample_particles( // nitadori's version
  const nbody_particle pb[], 
  const int nbody, 
  int sample_freq, 
  std::vector<float4> &sample_array, 
  int &nsample, 
  real &rmax)
{
  extern double get_ninteractions();
  double weight = get_ninteractions() / nbody * 1.e-4;
  if(0.0 == weight) weight = 1.0;
  sample_array.clear();
  for(int i=0,  ii=0; ii<nbody; i++, ii+=sample_freq){
    sample_array.push_back(float4(pb[ii].get_pos(), weight));
  }
  nsample = sample_array.size();
  MP_gather_sample_coords(sample_array);
  nsample = sample_array.size();

  if(0.0 == GLOBAL_ramx){
    double tx=0.0, ty=0.0, tz=0.0;
    for(int i=0; i<nbody; i++){
      vector3 r= pb[i].get_pos();
      tx = std::max(tx, fabs(r[0]));
      ty = std::max(ty, fabs(r[1]));
      tz = std::max(tz, fabs(r[2]));
    }
    rmax = std::max(tx, std::max(ty, tz));
    rmax = MP_doublemax(rmax);
    GLOBAL_ramx = rmax;
  }
  rmax = GLOBAL_ramx;
}

static int nmaxsample = NMAXSAMPLE;

static int determine_sample_freq(const int nbody)
{
    double nreal = nbody;
    MP_sum(nreal);
    double maxsample = (nmaxsample*0.8); // 0.8 is safety factor
    int sample_freq = int((nreal+maxsample)/maxsample);
    MP_int_bcast(sample_freq);
    return sample_freq;
}

static int sample_freq;
static std::vector<float4> sample_array;
static int nsample;
static int npx, npy, npz;

static inline int which_box(
  const vector3 &pos,
  const vector3 xlow[],
  const vector3 xhigh[]){
  int p = 0;
  if(pos[0] < xlow[p][0]) return -1;
  for(int ix=0; ix<npx; ix++, p+=npy*npz){
    if(pos[0] < xhigh[p][0]) break;
  }
  if(pos[0] > xhigh[p][0]) return -1;
  if(pos[1] < xlow[p][1]) return -1;
  for(int iy=0; iy<npy; iy++, p+=npz){
    if(pos[1] < xhigh[p][1]) break;
  }
  if(pos[1] > xhigh[p][1]) return -1;
  if(pos[2] < xlow[p][2]) return -1;
  for(int iz=0; iz<npz; iz++, p++){
    if(pos[2] < xhigh[p][2]) break;
  }
  if(pos[2] > xhigh[p][2]) return -1;
  return p;
}

int set_nmaxsample(int _nmaxsample){
  nmaxsample = _nmaxsample;
  MP_int_bcast(nmaxsample);
  sample_array.resize(nmaxsample);
  return 0;
}

void initialize_division(int nbody)
{
  static bool initcall = true;
  if(initcall){
    sample_freq = determine_sample_freq(nbody);
    create_division(MP_proccount(), npx, npy, npz);
    initcall = false;
  }
}

void setup_division(nbody_particle * pb, int nbody, int npdim[3], vector3* xlow, vector3* xhigh)
{
  static int initcall = true;
  if(initcall){
    initialize_division(nbody);
  }
  static real rmax;
  collect_sample_particles(pb, nbody, sample_freq, sample_array, nsample, rmax);
  rmax *= 1.2; // for skipping division
  npdim[0]=npx; npdim[1]=npy; npdim[2]=npz;
  if (MP_myprocid() == 0){
    determine_division_with_weight(nsample, sample_array, npx, npy, npz, rmax,xlow, xhigh);
  }
  int nwords=MP_proccount()*3;
  MP_double_bcast((real*)xlow, nwords);
  MP_double_bcast((real*)xhigh, nwords);
}

void exchange_particles_alltoall_vector(nbody_particle pb[],
  int &_nbody,
  const int nbmax,
  const vector3 xlow[],
  const vector3 xhigh[])
{
  const int nbody = _nbody;
  int myid = MP_myprocid();
  int nprocs = MP_proccount();
  static std::vector<nbody_particle> psend[NMAXPROC];
  static std::vector<nbody_particle> precv[NMAXPROC];
  bool initcall = true;
  if(initcall){
    initcall = false;
    for(int p=0; p<nprocs; p++){
      psend[p].reserve(64);
      precv[p].reserve(64);
    }
  }

  int iloc = 0;
  double dtsel = 1.e9;
  {
    RAII_Timer timer("Select particle (:p : ", &dtsel, std::cout, MP_root());
    const Boundary boundary(xlow[myid], xhigh[myid]);
    int i = 0;
    int j = nbody-1;
    for(;;){
      while( boundary.isinbox(pb[i].get_pos()) && i < nbody) i++;
      while(!boundary.isinbox(pb[j].get_pos()) && j >= 0) j--;
      if(i >= j) break;
      pb[i].swap(pb[j]);
      i++;
      j--;
    }
    iloc = i;

    for(int p=0; p<nprocs; p++){
      psend[p].clear();
      precv[p].clear();
    }
    for(int i=iloc; i<nbody; i++){
      int ibox = which_box(pb[i].get_pos(), xlow, xhigh);
      if(ibox < 0){
       cerr << myid <<" exchange_particle error: particle in no box..." << endl;
        vector3 fpos = pb[i].get_pos();
#ifdef PRINT
        unsigned long *upos = (unsigned long *)&fpos[0];
        cout << boost::format("[%f %f %f], [%lx %lx %lx]")
                % fpos[0] % fpos[1] % fpos[2]
                % upos[0] % upos[1] % upos[2]
             << endl;
#endif
        MP_Abort(1);
      }else{
        psend[ibox].push_back(pb[i]);
      }
    }
  }

  MP_sync();

  {
#ifdef PRINT
    double dtime = 1.e9;
    RAII_Timer timer("Exchange particle : ", &dtime, std::cout, MP_root());
#endif
    MP_alltoallv(psend, precv);
  }

  int nsendtot = 0, nrecvtot = 0;
  for(int p=0; p<nprocs; p++){
    nsendtot += psend[p].size();
    nrecvtot += precv[p].size();
  }
  const int nrecvloc = nrecvtot;
  double nbodytot = double(nbody);
  MP_sum(nbodytot);
  MP_int_sum(nsendtot);
  MP_int_sum(nrecvtot);
#ifdef PRINT
  double bw = 2.0 * double(sizeof(nbody_particle) * nsendtot) / dtime * 1.e-9;
  if(MP_root()){
    assert(nsendtot == nrecvtot);
    cout << "Exchanged particles = " << nsendtot 
         << ", " << nsendtot / nbodytot * 100 << "%"
         << ", " << dtime << "sec" << endl;
    cout << "Global Bandwidth " << bw << " GB/s" << endl;
  }
#endif
  if(!(iloc + nrecvloc <= nbmax)){
#ifdef PRINT
    printf("%d %d %d %d\n", myid, iloc, nrecvloc, nbmax);
#endif
    assert(iloc + nrecvloc <= nbmax);
  }
  for(int p=0; p<nprocs; p++){
    int size = precv[p].size();
    for(int i=0; i<size; i++){
      precv[p][i].copy(pb + iloc++);
    }
  }
  _nbody = iloc;
}

void exchange_particles_alltoall_vector_init(
  nbody_particle pb_src[],
  nbody_particle pb_dst[],
  int &nbody,
  const int nbmax,
  const vector3 xlow[],
  const vector3 xhigh[])
{
  int myid = MP_myprocid();
  int nprocs = MP_proccount();
  std::vector<nbody_particle> psend[NMAXPROC];
  std::vector<nbody_particle> precv[NMAXPROC];
  for(int p=0; p<nprocs; p++){
    psend[p].reserve(1024);
    precv[p].reserve(1024);
  }
  int idst = 0;
  const int NDIV = 64;
  for(int idiv=0; idiv<NDIV; idiv++){
#ifdef PRINT
    if(MP_root()){
      cout << "\nInitial exchange : " << idiv << "/" << NDIV << endl;
    }
#endif
    const int istart = ((0+idiv)*nbody)/NDIV;
    const int iend   = ((1+idiv)*nbody)/NDIV;

    for(int p=0; p<nprocs; p++){
      psend[p].clear();
      precv[p].clear();
    }
    for(int i=istart; i<iend; i++){
      int ibox = which_box(pb_src[i].get_pos(), xlow, xhigh);
      if(ibox < 0){
        cerr << myid <<" exchange_particle error: particle in no box..." << endl;
        vector3 fpos = pb_src[i].get_pos();
#ifdef PRINT
        unsigned long *upos = (unsigned long *)&fpos[0];
        cout << boost::format("[%f %f %f], [%lx %lx %lx]")
                % fpos[0] % fpos[1] % fpos[2]
                % upos[0] % upos[1] % upos[2]
             << endl;
#endif
        MP_Abort(1);
      }else{
        psend[ibox].push_back(pb_src[i]);
      }
    }

    {
#ifdef PRINT
      double dtime = 1.e9;
      RAII_Timer timer("Exchange particle : ", &dtime, std::cout, MP_root());
#endif
      MP_alltoallv(psend, precv);
    }

    int nsendtot = 0, nrecvtot = 0;
    for(int p=0; p<nprocs; p++){
      nsendtot += psend[p].size();
      nrecvtot += precv[p].size();
    }
    MP_int_sum(nsendtot);
    MP_int_sum(nrecvtot);
#ifdef PRINT
    double bw = 2.0 * double(sizeof(nbody_particle) * nsendtot) / dtime * 1.e-9;
    if(MP_root()){
      assert(nsendtot == nrecvtot);
      cout << "Exchanged particles = " << nsendtot << ", " << dtime << "sec" << endl;
      cout << "Global Bandwidth " << bw << " GB/s" << endl;
    }
#endif
    for(int p=0; p<nprocs; p++){
      int size = precv[p].size();
      for(int i=0; i<size; i++){
        precv[p][i].copy(pb_dst + idst++);
        assert(idst <= nbmax);
      }
    }
  } // idiv
  nbody = idst;
  {
    Boundary boundary(xlow[myid], xhigh[myid]);
    for(int i=0; i<nbody; i++){
      assert(boundary.isinbox(pb_dst[i].get_pos()));
    }
  }
}
