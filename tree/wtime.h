#ifndef WTIME_H
#define WTIME_H
#include <sys/time.h>

static double get_wtime(){
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + 1.e-6 * tv.tv_usec;
}

#include <string>
#include <iostream>

struct RAII_Timer{
  const std::string name;
  double * const tout;
  std::ostream &os;
  const bool outflag;
  double t0;
  RAII_Timer(
    const std::string &_name, 
    double            *_tout    = NULL, 
    std::ostream      &_os      = std::cout, 
    bool               _outflag = true) :
    name(_name), tout(_tout), os(_os), outflag(_outflag), t0(get_wtime())
  {
  }
  ~RAII_Timer(){
    double t1 = get_wtime();
#ifdef PRINT
    if(outflag) os << name << t1-t0 << " sec" << std::endl;
#endif
    if(tout) *tout = t1-t0;
  }
};
#endif
