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

/*! \file gpu_generic_integrator 
 *  \brief kernel code for gpu_generic_integrator
 *
*/

#include "swarm.h"
#include <cux/cux.h>
#include <astro/util.h>
#include <cassert>
#include "swarmlog.h"

namespace swarm {

/*!
	\brief gpu_generic_integrator - Versatile GPU integrator template

	gpu_generic_integrator is a template class designed to make it easy to implement
	powerful GPU-based integrators by supplying a 'propagator' class (that provides
	a function to advance the ensemble by one timestep), and a 'stopper' class (that
	tests whether the integration should stop for a particular system in the
	ensemble). Given the two classes, gpu_generic_integrator acts as a driver
	routine repeatedly calling the GPU kernel until all systems in the ensemble have
	finished the integration. It does this in an optimized manner, taking care to
	compactify() the ensemble when needed to efficiently utilize GPU resources
	(NOTE: compactification not yet implemented).

	For a canonical example of how to build an integrator using
	gpu_generic_integrator, look at the gpu_euler integrator.

	Integration loop outline (I == gpu_generic_integrator object, H == propagator
	object, stop == stopper object):

	I.integrate(ens):
		if(ens.last_integrator() != this):
			H.initialize(ens);
			stop.initialize(ens);
		do:
			gpu_integrate<<<>>>(max_step, ens, (gpu_t)H, (gpu_t)stop)	[m'threaded execution on the GPU]
				while step < max_step:
					if(T < Tend || stop()):
						ens(sys).flags |= INACTIVE
						break;
					H.advance():			[ implementation supplied by the developer ]
						foreach(bod in sys):
							advance bod
							call stop.test_body
						return new time T
			if(beneficial):
				ens.compactify();
		while(active_systems > 0)

	To build an integrator using gpu_generic_integrator, the developer must supply a
	propagator and a stopper class that conform to the following interfaces:

	// propagator class: advance the system by one time step
	//
	// CPU state and interface. Will be instantiated on construction
	// of gpu_generic_integrator object. Keep any data that need
	// to reside on the CPU here.
	struct propagator
	{
		// GPU state and interface (per-grid). Will be passed 
		// as an argument to integration kernel. Any per-block read-only
		// variables should be members of this structure.
		struct gpu_t
		{
			// GPU per-thread state and interface. Will be instantiated
			// in the integration kernel. Any per-thread
			// variables should be members of this structure.
			struct thread_state_t
			{
				__device__ thread_state_t(const gpu_t &H, ensemble &ens, const int sys, double T, double Tend);
			};

			// Advance the system - this function must advance the system
			// sys by one timestep, making sure that T does not exceed Tend.
			// Must return the new time of the system.
			//
			// This function MUST also call stop.test_body() for every body
			// in the system, after that body has been advanced by a timestep.
			template<typename stop_t>
			__device__ double advance(ensemble &ens, thread_state_t &pt, int sys, double T, double Tend, stop_t &stop, typename stop_t::thread_state_t &stop_ts, int step)
		};

		// Constructor will be passed the cfg object with the contents of
		// integrator configuration file. It will be called during construction
		// of gpu_generic_integrator. It should load any values necessary
		// for initialization.
		propagator(const config &cfg);

		// Initialize temporary variables for ensemble ens. This function
		// should initialize any temporary state that is needed for integration
		// of ens. It will be called from gpu_generic_integrator, but only
		// if ens.last_integrator() != this. If any temporary state exists from
		// previous invocation of this function, it should be deallocated and
		// the new state (re)allocated.
		void initialize(ensemble &ens);

		// Cast operator for gpu_t. This operator must return the gpu_t object
		// to be passed to integration kernel. It is called once per kernel
		// invocation.
		operator gpu_t();
	};


	// stopper class: mark a system inactive if conditions are met
	//
	// CPU state and interface. Will be instantiated on construction
	// of gpu_generic_integrator object. Keep any data that need
	// to reside on the CPU here.
	struct stopper
	{
		// GPU state and interface (per-grid). Will be passed 
		// as an argument to integration kernel. Any per-block read-only
		// variables should be members of this structure.
		struct gpu_t
		{
			// GPU per-thread state and interface. Will be instantiated
			// in the integration kernel. Any per-thread
			// variables should be members of this structure.
			struct thread_state_t
			{
				__device__ thread_state_t(gpu_t &stop, ensemble &ens, const int sys, double T, double Tend);
			};

			// test any per-body stopping criteria for body (sys,bod). If 
			// your stopping criterion only depends on (x,v), test for it 
			// here. This will save you the unnecessary memory accesses 
			// that would otherwise be made if the test was made from 
			// operator().
			//
			// Called _after_ the body 'bod' has advanced a timestep.
			//
			// Note: you must internally store the result of your test,
			// and use/return it in subsequent call to operator().
			//
			__device__ void test_body(thread_state_t &ts, ensemble &ens, int sys, int bod, double T, double x, double y, double z, double vx, double vy, double vz);

			// Called after a system sys has been advanced by a timestep.
			// Must return true if the system sys is to be flagged as
			// INACTIVE (thus stopping further integration)
			__device__ bool operator ()(thread_state_t &ts, ensemble &ens, int sys, int step, double T);
		};

		// Constructor will be passed the cfg object with the contents of
		// integrator configuration file. It will be called during construction
		// of gpu_generic_integrator. It should load any values necessary
		// for initialization.
		stopper(const config &cfg);

		// Initialize temporary variables for ensemble ens. This function
		// should initialize any temporary state that is needed for integration
		// of ens. It will be called from gpu_generic_integrator, but only
		// if ens.last_integrator() != this. If any temporary state exists from
		// previous invocation of this function, it should be deallocated and
		// the new state (re)allocated.
		void initialize(ensemble &ens);

		// Cast operator for gpu_t. This operator must return the gpu_t object
		// to be passed to integration kernel. It is called once per kernel
		// invocation.
		operator gpu_t();
	};

*/

/********* CLASS DEFINITION *********/
template<typename stopper_t, typename propagator_t>
class gpu_generic_integrator : public integrator
{
protected:
	stopper_t	stop;
	propagator_t	H;

	int gpu_ensemble_id;
	int steps_per_kernel_run;

	dim3 gridDim;
	int threadsPerBlock;

	/// temp variable for return values (gpu pointer)
	cuxDeviceAutoPtr<retval_t> retval_gpu;	

public:
	gpu_generic_integrator(const config &cfg);

public:
	void integrate(gpu_ensemble &ens, double T);
};

/*!
 \brief gpu generic integrator 
 ***************** CONSTRUCTOR*****************

 ...
 @param[in] cfg
*/
template<typename stopper_t, typename propagator_t>
gpu_generic_integrator<stopper_t, propagator_t>::gpu_generic_integrator(const config &cfg)
        : H(cfg), stop(cfg), retval_gpu(1)
{
        steps_per_kernel_run = cfg.count("steps per kernel run") ? static_cast<int>(std::floor(atof(cfg.at("steps per kernel run").c_str()))) : 1000;
        threadsPerBlock = cfg.count("threads per block") ? atoi(cfg.at("threads per block").c_str()) : 64;
        gpu_ensemble_id = cfg.count("gpu_ensemble_id") ? atoi(cfg.at("gpu_ensemble_id").c_str()) : 0;

        // test the number of threads for validity
        if(threadsPerBlock <= 0) { ERROR("'threads per block' must be greater than zero (currently, " + str(threadsPerBlock) + ")"); }
        if(!is_power_of_two(threadsPerBlock)) { ERROR("'threads per block' must be a power of two (currently, " + str(threadsPerBlock) + ")"); }
        if(threadsPerBlock > MAXTHREADSPERBLOCK) { ERROR("'threads per block' cannot be greater than " + str(MAXTHREADSPERBLOCK) + " (currently, " + str(threadsPerBlock) + ")"); }
}

/********************* define integrate function here **************
 *
 * integrate initializes the propagator and stopper and then executes 
 * the specified kernel using the gpu_integ_driver function in swarmlib.cu
 *
 */


/*!
 \brief gpu integrate 
 
 @param[out] ens
 @param[in] dT
*/
template<typename stopper_t, typename propagator_t>
void gpu_generic_integrator<stopper_t, propagator_t>::integrate(gpu_ensemble &ens, double dT)
{
	// Upload the kernel parameters
	if(ens.last_integrator() != this)
	{
		ens.set_last_integrator(this);
		configure_grid(gridDim, threadsPerBlock, ens.nsys());

		// upload ensemble
		assert((gpu_ensemble_id>=0)&&(gpu_ensemble_id<MAX_GPU_ENSEMBLES));
		cuxErrCheck( cudaMemcpyToSymbol(gpu_integ_ens[gpu_ensemble_id], &ens, sizeof(gpu_integ_ens[gpu_ensemble_id])) );

		// initialize propagator, stopping condition
		H.initialize(ens);
		stop.initialize(ens);

		if(dT == 0.) { return; }
	}

	// flush CPU/GPU logs
	log::flush(log::memory | log::if_full);

	// initialize stop-time temporary array
	cuxDeviceAutoPtr<real_time> Tstop(ens.nsys());

	// execute the kernel in blocks of steps_per_kernel_run timesteps
	int nactive0 = -1;
	int iter = 0;
	do
	{
//		lprintf(hlog, "Starting kernel run #%d\n", iter);
//		lprintf(hlog, "Another unnecessary message from the CPU side\n");

		retval_gpu.memset(0);
		gpu_integ_driver<typename stopper_t::gpu_t, typename propagator_t::gpu_t><<<gridDim, threadsPerBlock>>>(retval_gpu, steps_per_kernel_run, H, stop, gpu_ensemble_id, Tstop, iter == 0 ? dT : 0.);
		cuxErrCheck( cudaThreadSynchronize() );
		iter++;

		// flush CPU/GPU logs
		log::flush(log::memory | log::if_full);

		retval_t retval;
		retval_gpu.get(&retval);
		if(nactive0 == -1) { nactive0 = retval.nrunning; }

		// check if we should compactify or stop
		if(retval.nrunning == 0)
		{
			break;
		}
		if(retval.nrunning - nactive0 > 128)
		{
			// TODO: compactify here
			nactive0 = retval.nrunning;
		}
	} while(true);
//	lprintf(hlog, "Exiting integrate");

	log::flush(log::memory);
}

} // end namespace swarm
