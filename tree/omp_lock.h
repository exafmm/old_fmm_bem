// Copyright(C) 2010 by 
// Tsuyoshi Hamada <hamada@progrape.jp>
// Keigo Nitadori <nitadori@margaux.astron.s.u-tokyo.ac.jp>
// Rio Yokota <rio.yokota@bristol.ac.uk>

#include <omp.h>
#ifdef LOCK_GPU
struct omp_lock{
	omp_lock_t lock;
	omp_lock(){ omp_init_lock(&lock); }
	~omp_lock() { omp_destroy_lock(&lock); }
	void set_lock() { omp_set_lock(&lock); }
	void unset_lock() { omp_unset_lock(&lock); }
	int test_lock() { return omp_test_lock(&lock); }
};
#else
struct omp_lock{
	omp_lock(){}
	~omp_lock() {}
	void set_lock() {}
	void unset_lock() { }
	int test_lock() { return 0; }
};
#endif
