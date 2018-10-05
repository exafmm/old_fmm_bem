/*
 *  time.C: Interface to some standard UNIX time functions.
 *
 *          WARNING! -- the "standard" functions are anything but standard,.
 *                      They may be both non-portable and subject to change!
 *
 *                      The code below is tested only in SunOS 4.1.3 and HPUX.
 *
 *.............................................................................
 *    version 1:  Mar 1994   Steve McMillan
 *    version 2:  
 *.............................................................................
 *  non-local functions: 
 *    cpu_time
 *.............................................................................
 */
#define real double
#include "stdlib.h"

#include <sys/times.h>
    struct tms buffer;
#include <unistd.h>
#include <sys/timeb.h>
#include <sys/timeb.h>
#include <sys/time.h>
   struct timeb wall_buffer;
extern "C" int ftime(struct timeb *tp);

  static long ticks_per_sec = 0;
  static real initial_cpu = 0;

//-----------------------------------------------------------------------------
//  cpu_init: initialize the CPU timer.
//-----------------------------------------------------------------------------

void cpu_init()
{
    times(&buffer);

    // Use both system and user time because of ambiguities
    // with Linux multiprocessing...

    initial_cpu = (real) (buffer.tms_utime + buffer.tms_stime);

    ticks_per_sec = sysconf(_SC_CLK_TCK);	// Clock ticks per second
}

//-----------------------------------------------------------------------------
//  cpu_time: return total processor time (in seconds) used since the timer
//            was last initialized by a call to cpu_init.
//-----------------------------------------------------------------------------

real cpu_time()
{

    if (!ticks_per_sec)
	cpu_init();

    times(&buffer);
    return ((real) (buffer.tms_utime + buffer.tms_stime - initial_cpu))
      			/ ticks_per_sec;
}


//-----------------------------------------------------------------------------
//  wall_init: initialize the wallclock timer.
//-----------------------------------------------------------------------------

static real initial_wall = 0;
void wall_init()
{
    ftime(&wall_buffer);

    initial_wall = ((real) wall_buffer.time)
	+((real) wall_buffer.millitm)/1000;
}
real wall_time()
{
    ftime(&wall_buffer);

    return  ((real) wall_buffer.time)
	+((real) wall_buffer.millitm)/1000 - initial_wall;
}

