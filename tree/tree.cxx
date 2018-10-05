// Copyright(C) 2010 by 
// Tsuyoshi Hamada <hamada@progrape.jp>
// Keigo Nitadori <nitadori@margaux.astron.s.u-tokyo.ac.jp>
// Rio Yokota <rio.yokota@bristol.ac.uk>

#include "../include/constants.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <cassert>
#include <fstream>
#include <iostream>
#include <omp.h>
using namespace std;

#define real double
#include "BHtree.h"
#include "nbody.h"
#include "nbody_particle.h"
#include "vector3.h"
#include "wtime.h"

#define PR(x)  cerr << #x << " = " << x << " "
#define PRC(x) cerr << #x << " = " << x << ",  "
#define PRL(x) cerr << #x << " = " << x << "\n"

#define MORTON_SKIP   1
#define DIVISION_SKIP 1

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi,*vi;
extern float *xj,*yj,*zj,*gxj,*gyj,*gzj,*vj;

typedef nbody_particle real_particle;
typedef nbody_system real_system;
typedef nbody_VF_ptr real_VF_ptr;
typedef nbody_RF_ptr real_RF_ptr;
typedef nbody_RRF_ptr real_RRF_ptr;

void real_system::set_nsize(const int _nsize){
		nsize = _nsize;
		pbuf1.resize(nsize);
		pbuf2.resize(nsize);
		log_allocate("real particle", nsize * sizeof(real_particle));
		log_allocate("real particle", nsize * sizeof(real_particle));
}

void real_system::gen_uniform(int nuniform){
	this->n = nuniform;
	this->time = 0.0;
	int nsize = int(1.25 * this->n);
	this->set_nsize(nsize);
        pb = this->get_pbuf1();

	for(int i=0; i<this->n; ){
		pb[i].set_index(long(this->n) * MP_myprocid() + i);
		pb[i].set_mass(vj[i]);
		pb[i].set_pos(vector3(xj[i], yj[i], zj[i]));
		pb[i].set_vel(vector3(0.0, 0.0, 0.0));
		i++;
	}
}

const real piinv = 1.0/M_PI;

void real_system::apply_vf(real_VF_ptr f)
{
#pragma omp parallel for
	for(int i = 0;i<n;i++){
		(pb[i].*f)();
	}
}

void tree(int& ni, int& nj, int neq)
{
        static real_system pb;
        int nmaxsample = NMAXSAMPLE;
        real theta = 0.4;
        int ncrit = 2000;
        int nuniform = ni;
        wall_init();
        set_nmaxsample(nmaxsample);
        MP_int_bcast(nuniform);
        initialize_division(pb.n);
        pb.gen_uniform(nuniform);
        assert(pb.n);
        static vector3 *xlow  = pb.xlowp;
        static vector3 *xhigh = pb.xhighp;
        vector3 r,acc;
        int npdim[3];
        setup_division(pb.get_particle_pointer(),pb.n,npdim,xlow,xhigh);
        MP_sync();
        if(MP_proccount() > 1){
                exchange_particles_alltoall_vector_init(
                        pb.get_pbuf1(), pb.get_pbuf2(),
                        pb.n, pb.get_nsize(), xlow, xhigh);
                pb.swap_pbuf();
        }
        setup_division(pb.get_particle_pointer(),pb.n,npdim,xlow,xhigh);
        exchange_particles_alltoall_vector(pb.get_particle_pointer(),pb.n,pb.get_nsize(), xlow, xhigh);
        MP_sync();
        pb.eps2_for_gravity = eps;
        pb.theta_for_tree = theta;
        pb.ncrit_for_tree = ncrit;
        pb.do_morton_sort = true;
        MP_sync();
        pb.calculate_gravity();
        for(int i=0; i<pb.n; i++) {
          if(MP_proccount() > 1){
            r = pb.get_pbuf1()[i].get_pos();
            acc = pb.get_pbuf1()[i].get_acc_gravity();
            vj[i] = pb.get_pbuf1()[i].get_mass();
            vi[i] = -0.25/pi*pb.get_pbuf1()[i].get_phi_gravity();
          } else {
            r = pb.get_pbuf2()[i].get_pos();
            acc = pb.get_pbuf2()[i].get_acc_gravity();
            vj[i] = pb.get_pbuf2()[i].get_mass();
            vi[i] = -0.25/pi*pb.get_pbuf2()[i].get_phi_gravity();
          }
          xj[i] = r[0];
          yj[i] = r[1];
          zj[i] = r[2];
          xi[i] = xj[i];
          yi[i] = yj[i];
          zi[i] = zj[i];
          gxi[i] = 0.25/pi*acc[0];
          gyi[i] = 0.25/pi*acc[1];
          gzi[i] = 0.25/pi*acc[2];
        }
        ni = pb.n;
        nj = pb.n;
#pragma omp parallel
	{
		vforce_close(omp_get_thread_num());
	}
}
