/*************************************************************************
 * Copyright (C) 2008-2010 by Mario Juric & Swarm-NG Development Team    *
 *                                                                       *
 * This program is free software; you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation; either version 3 of the License.        *
 *                                                                       *
 * This program is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 * GNU General Public License for more details.                          *
 *                                                                       *
 * You should have received a copy of the GNU General Public License     *
 * along with this program; if not, write to the                         *
 * Free Software Foundation, Inc.,                                       *
 * 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ************************************************************************/

/*! \file swarmlib.cu
 *  \brief kernel code for gpu_generic_integrator
 *
*/

#include "swarm.h"
#include <cux/cux.h>
#include <astro/util.h>
#include <cassert>
#include "swarmlog.h"
#include "ThreeVector.hpp"


gpulog::host_log hlog;

#if __CUDACC__
// The assumption is all CUDA code will be concatenated/included and compiled
// as a single source file (thus avoiding the creation of duplicate copies of 
// hlog and dlog)
__constant__ gpulog::device_log dlog;
#endif


/*********** USEFUL UTILITY FUNCTIONS ****************/
#if __DEVICE_EMULATION__
	extern "C" void debug_hook();
#else
	#define debug_hook()
#endif

// ------------------------ THREAD functions -------------------------------------- //

/*!
  \brief Computes the global linear ID of the thread. Used from kernels.

  NOTE: Supports 3D grids with 1D blocks of threads
  @return threadId
*/
inline __device__ uint32_t threadId()
{
// This will be in inner loops, so may want to optimize
#if USE_1D_GRID
	const uint32_t id = blockIdx.x * blockDim.x + threadIdx.x;
#else
	const uint32_t id = ((blockIdx.z * gridDim.y + blockIdx.y) * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;
#endif
	return id;
}

/*!
  \brief Computes the number of threads per block for 1D or 3D grid. Used from kernels.

  NOTE: Supports 3D grids with 1D blocks of threads
  @return threadsPerBlock
*/
inline __device__ uint32_t threadsPerBlock()
{
#if USE_1D_GRID
	return blockDim.x;
#else
	return blockDim.x * blockDim.y * blockDim.z;
#endif
}


namespace swarm {

/************************* UTILITY FUNCTIONS **************************/
///////////////
//// Some generically useful force/jerk computation functions
///////////////

/*!
  Calculate accelerations for the system, storing the output into array aa

   NOTE: The second loop goes from (nbod..0], which is the optimal choice
   from numerical precision standpoint if bod=0 is the most massive
   in the system.

  @param[in] ens
  @param[in] sys
  @param[out] aa
*/
template<typename real>
__device__ void compute_acc(ensemble &ens, const int sys, const cuxDevicePtr<real, 3> &aa)
{
        typedef ThreeVector<real> V3;

        for ( unsigned int i=0;i<ens.nbod();++i )
        {
                V3 xi( ens.x ( sys,i ), ens.y ( sys,i ), ens.z ( sys,i ) );
                V3 ai(0.);
                for (int j=ens.nbod()-1; j >= 0; j--)
                {
                        if ( j==i ) continue; // Ignore body interacting with itself

                        V3 dx(ens.x(sys,j), ens.y(sys,j), ens.z(sys,j));  dx -= xi;
                        real r2 = dx.MagnitudeSquared();
                        real rinv = rsqrt ( r2 );
                        rinv *= ens.mass ( sys,j );
                        real rinv3 = rinv/r2;

                        dx *= rinv3;
                        ai += dx;
                }
                aa ( sys, i, 0 ) = ai.X();
                aa ( sys, i, 1 ) = ai.Y();
                aa ( sys, i, 2 ) = ai.Z();
        } // end loop over bodies
}

/*!
  Calculate acceleration and jerk for the system, storing the outputs into arrays aa, jj
  Is the ordering of dimensions what we want?

        NOTE: The second loop goes from (nbod..0], which is the optimal choice
        from numerical precision standpoint if bod=0 is the most massive
        in the system.

  @param[in] ens
  @param[in] sys
  @param[out] aa
  @param[out] jj
*/
template<typename real>
__device__ void compute_acc_jerk(ensemble &ens, const int sys, const cuxDevicePtr<real, 3> &aa, const cuxDevicePtr<real, 3> &jj)
{
        typedef ThreeVector<real> V3;

        for ( unsigned int i=0;i<ens.nbod();++i )
        {
                V3 xi( ens.x ( sys,i ), ens.y ( sys,i ), ens.z ( sys,i ) );
                V3 vi( ens.vx ( sys,i ),ens.vy ( sys,i ),ens.vz ( sys,i ) );
                V3 ai(0.), ji(0.);
                for (int j=ens.nbod()-1; j >= 0; j--)
                {
                        if ( j==i ) continue; // Ignore body interacting with itself
                        V3 dx(ens.x(sys,j), ens.y(sys,j), ens.z(sys,j));  dx -= xi;
                        V3 dv(ens.vx(sys,j),ens.vy(sys,j),ens.vz(sys,j)); dv -= vi;
                        real r2 = dx.MagnitudeSquared();
                        real rv = dot ( dx,dv );
                        real rinv = rsqrt ( r2 );
                        rv *= 3./r2;
                        rinv *= ens.mass ( sys,j );
                        real rinv3 = rinv/r2;

                        dx *= rinv3;
                        ai += dx;
                        dv *= rinv3;
                        ji += dv;
                        dx *= rv;
                        ji -= dx;
                }

                aa ( sys, i, 0 ) = ai.X();
                aa ( sys, i, 1 ) = ai.Y();
                aa ( sys, i, 2 ) = ai.Z();
                jj ( sys, i, 0 ) = ji.X();
                jj ( sys, i, 1 ) = ji.Y();
                jj ( sys, i, 2 ) = ji.Z();
        } // end loop over bodies
}

//---------------------------------------------------------------

/*!
  \brief empty kernel to simplify setting of L1 cache size
*/
__global__ void dummy() {};

/*!
  \brief Sets larger avaliable L1 cache size to be default (48k on GF100)
*/
void set_cuda_cache_large()
{
cudaFuncSetCacheConfig(dummy, cudaFuncCachePreferL1);
dummy<<<16,32>>>();
};

/*!
  \brief Sets smaller avaliable L1 cache size to be default (16k on GF100)
*/
void set_cuda_cache_small()
{
cudaFuncSetCacheConfig(dummy, cudaFuncCachePreferShared);
dummy<<<16,32>>>();
};

/*!
  \brief compute and store the number of systems that remain active

  NOTE: assumes not more than MAXTHREADSPERBLOCK threads per block
  NOTE: assumes a nthreads is a power of 2
  NOTE: assumes *nrunning = 0 on input
 @param[out] nrunning number of active systems
 @param[in] ens ensemble
*/
static const int MAXTHREADSPERBLOCK = 256;
__device__ void count_running(int *nrunning, double *Tstop, ensemble &ens)
{
#if 0   // If you don't want to use shared memory for this
	// Direct counting (via atomicAdd) of the number of running threads. Slower
	// but requires no shared memory.
	// TODO: Test if this implementation is fast enough to be activated
	int sys = threadId();
	int running = sys < ens.nsys() ? !(ens.flags(sys) & ensemble::INACTIVE || ens.time(sys) >= Tstop[sys]) : 0;
	if(running) { atomicAdd(nrunning, 1); /*printf("sys=%d running.\n", sys);*/ }
	return;
#else
	// We think it's ok to make this extern ?
	__shared__ int running[MAXTHREADSPERBLOCK];	// takes up 1k of shared memory (for MAXTHREADSPERBLOCK=256)

	int tpb = threadsPerBlock();
	int sys = threadId();
	const int widx = sys % tpb;
	running[widx] = sys < ens.nsys() ? !(ens.flags(sys) & ensemble::INACTIVE || ens.time(sys) >= Tstop[sys]) : 0;

	// Prefix sum algorithm (assumes block size <= MAXTHREADSPERBLOCK).
	// 1) sum up the number of running threads in this block
	for(int i = 2; i <= tpb; i *= 2)
	{
		 __syncthreads();
		if(widx % i == 0) { running[widx] += running[widx + i/2]; }
	}

	// 2) add the sum to the total number of running threads
	if(widx == 0)
	{
		atomicAdd(nrunning, running[0]);
	}
#endif
}



#if 0
/*!
  \brief Standalone kernel for counting the number of active ensembles in a GPU ensemble

 @param[out] nrunning number of active systems 
 @param[in] ens ensemble
*/
__global__ void get_nactive_kernel(int *nrunning, ensemble ens)
{
	count_running(nrunning, ens);
}
#endif

/// return true if input parameter is a power of two
inline bool is_power_of_two(int x) { return !(x & (x-1)); }


/***************************** MAIN BODY OF GPU INTEGRATOR CODE **********************
 *
 *   gpu_integ_driver is called from gpu_generic_integrator::integrate.  It in turn
 *   calls generic_integrate_system.
 *
 *   generic_integrate_system is a generic versatile framework that will implement
 *   the methods in the specific propagator and stopper to do the actual work.  The
 *   calculations are actually done in each propagator's advance method.
 *
 */
 
///////////////////
//// Generic versatile integrator template framework
///////////////////

#define MAX_GPU_ENSEMBLES 4
__constant__ ensemble gpu_integ_ens[MAX_GPU_ENSEMBLES];

/// generic GPU integrator return value
struct retval_t
{
	int nrunning;
};

/*!
 \brief generic integrate system  
 
 ...
 ...
 @param[out] retval
 @param[in,out] ens
 @param[in] sys
 @param[in] max_steps
 @param[in,out] H
 @param[in,out] stop
*/
template<typename stopper_t, typename propagator_t>
__device__ void generic_integrate_system(retval_t *retval, ensemble &ens, int sys, int thr, int max_steps, propagator_t &H, stopper_t &stop, double Tstop)
{
	// initialize propagator and stopper per-thread states
	double T = ens.time(sys);
	typename stopper_t::thread_state_t    stop_ts(stop, ens, sys, T, Tstop);
	typename propagator_t::thread_state_t    H_ts(H, ens, sys, T, Tstop);

	// advance the system until we reach max_steps, Tstop, or stop becomes true
	unsigned int step = 0;
	while(true)
	{
		// stopping conditions
		if(stop_ts.is_step_complete())
  		   {
		   if(T >= Tstop) 				{ if(T >= ens.time_end(sys)) { ens.flags(sys) |= ensemble::INACTIVE; } break; }
		   if(stop(stop_ts, ens, sys, step, T)) 	{ ens.flags(sys) |= ensemble::INACTIVE; break; }
		   if(step == max_steps) 			{ break; }

  		   log::output_system_if_needed(dlog, ens, T, sys);
 	 	   }

		// actual work
		T = H.advance(ens, H_ts, sys, thr, T, Tstop, stop, stop_ts, step);

		if(stop_ts.is_step_complete())
   		   step++;
	}
	ens.nstep(sys) += step;
	log::output_system_if_needed(dlog, ens, T, sys);

	ens.time(sys) = T;
}

/*!
 \brief gpu integrate driver 
 
 ...
 @param[out] retval
 @param[in] max_steps
 @param[in] H
 @param[in] stop
 @param[in] gpu_ensemble_id
*/
template<typename stopper_t, typename propagator_t>
__global__ void gpu_integ_driver(retval_t *retval, int max_steps, propagator_t H, stopper_t stop, const int gpu_ensemble_id, real_time *Tstop, double dT)
{
	// Need to test that don't get bogus gpu_ensemble_id
 	if((gpu_ensemble_id<0)||(gpu_ensemble_id>=MAX_GPU_ENSEMBLES))
		// Is this proper way to bail from GPU?
 		{ retval->nrunning = -1; return; }	

	// find the system we're to work on
	ensemble &ens = gpu_integ_ens[gpu_ensemble_id];
	int tps = propagator_t::threads_per_system(ens.nbod());
	int sys = threadId() / tps;
	int thr = threadId() % tps;

#if __DEVICE_EMULATION__
	if(sys == ens.nsys()-1)
	{
		double pct = 100. * (double)dlog.size() / dlog.capacity();
		printf("sys=%d T=%g Tend=%g Tstop=%g dT=%g flags=%d max_steps=%d logpct=%.2g%% [%d/%d]\n", sys, ens.time(sys), ens.time_end(sys), Tstop[sys], dT, ens.flags(sys), max_steps, pct, dlog.size(), dlog.capacity());
	}
#endif
	if(sys < ens.nsys() && !(ens.flags(sys) & ensemble::INACTIVE))
	{
		if(dT != 0.)
		{
			// if dT != 0, recompute the stop time in Tstop array as the current 
			// time + dT, making sure not to exceed the system end time, ens.time_end().

			Tstop[sys] = min(ens.time_end(sys), ens.time(sys) + dT);
		}

		generic_integrate_system(retval, ens, sys, thr, max_steps, H, stop, Tstop[sys]);
	}

#if __DEVICE_EMULATION__
	if(sys == ens.nsys()-1)
	{
		double pct = 100. * (double)dlog.size() / dlog.capacity();
		printf("sys=%d T=%g Tend=%g Tstop=%g dT=%g flags=%d max_steps=%d logpct=%8.4g%% [%d/%d]\n", sys, ens.time(sys), ens.time_end(sys), Tstop[sys], dT, ens.flags(sys), max_steps, pct, dlog.size(), dlog.capacity());
	}
#endif

	count_running(&retval->nrunning, Tstop, ens);
	
#if __DEVICE_EMULATION__
	if(sys == ens.nsys()-1)
	{
		printf("sys=%d nrunning=%d\n", sys, retval->nrunning);
	}
#endif
}

// Templated acc_jerk calculation functions


} // end namespace swarm
