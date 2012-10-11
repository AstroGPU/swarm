/*************************************************************************
 * Copyright (C) 2011 by Saleh Dindar and the Swarm-NG Development Team  *
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

/*! \file bppt.hpp
 *   \brief Defines the GPU integrator class with one thread for each body-pair.
 *          also defines the interface to GPU and CUDA. 
 *
 *
 */


#pragma once

#include "../integrator.hpp"
#include "../log/log.hpp"
#include "../plugin.hpp"

#include "helpers.hpp"
#include "utilities.hpp"
#include "device_settings.hpp"


namespace swarm {
namespace gpu {

/*! Class of GPU integrators with a thread for each body-pair
 * \addtogroup integrators
 *
 *  Using a thread for each body-pair is to parallelize as much as possible
 *  when integrating an ensemble. The thread assignment is as follows
 *  1. When computing interaction forces (and higher derivatives) between bodies, one thread is
 *  assigned to each pair of bodies.
 *  2. When integrating quantities for bodies individually, one thread is assigned to each 
 *  coordinate component of each body.
 *  3. When advancing the time or checking for stop criteria or setting the time step, only
 *  one thread is used.
 *
 *  We use predicate barriers for each case since the number of threads that are actually working
 *  is not the same in each case. For example, in case of 3 bodies, there are 3 pairs, 9 body-components and
 *  1 thread is needed for overall tasks.
 *
 *  For better coalesced reads from global and shared memory, the block is structured in a 
 *  non-traditional way. Innermost (x) component is part of the system id, and the other component
 *  is the body-pair.
 *
 *  Global functions defined here are used inside kernels for consistent interpretation 
 *  of thread and shared memory references
 *
 *  The integrators that derive from this may override three functions:
 *  thread_per_system(), shmem_per_sytem(), 
 *
 */
namespace bppt {


/**
 * Kernel Helper Function: Extract system ID from CUDA thread ID
 */
GPUAPI int sysid(){
	return ((blockIdx.z * gridDim.y + blockIdx.y) * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;
}

/**
 * Kernel Helper Function: Extract system sequence number inside current block
 */
GPUAPI int sysid_in_block(){
	return threadIdx.x;
}

/**
 * Kernel Helper Function: Extract the worker-thread number for current system
 */
GPUAPI int thread_in_system() {
	return threadIdx.y;
}


/**
 * Kernel Helper Function: Extract number of systems per a block from CUDA thread information.
 */
GPUAPI int system_per_block_gpu() {
    return blockDim.x;
  };

/**
 * Kernel Helper Function: Logical coordinate component id [1:x,2:y,3:z] calculated from thread ID info
 */
GPUAPI int thread_component_idx(int nbod) {
	return thread_in_system() / nbod;
}

/**
 * Kernel Helper Function: Logical body id [0..nbod-1] calculated from thread ID info
 */
GPUAPI int thread_body_idx(int nbod) {
	return thread_in_system() % nbod;
}


/**
 * Kernel Helper Function: Get the pointer to dynamic shared memory allocated for the system.
 * This function assumes that the memory is used through CoalescedStructArray with a chunk size
 * of SHMEM_CHUNK_SIZE. This uses overlapping data structures to provide coalescing for shared
 * memory.
 */
template< class Impl, class T> 
GPUAPI void * system_shared_data_pointer(Impl* integ, T compile_time_param) {
	extern __shared__ char shared_mem[];
	int b = sysid_in_block() / SHMEM_CHUNK_SIZE ;
	int i = sysid_in_block() % SHMEM_CHUNK_SIZE ;
	int idx = i * sizeof(double) 
		+ b * SHMEM_CHUNK_SIZE 
		* Impl::shmem_per_system(compile_time_param);
	return &shared_mem[idx];
}



/**
 *  \brief Common functionality and skeleton for body-pair-per-thread integrators
 *  Common tasks include:
 *   - Setting up the CUDA grid and block values.
 *   - Calculating the amount of shared memory.
 */
class integrator : public gpu::integrator  {
	typedef gpu::integrator Base;
	protected:
	//! Number of systems allocated in a block. 
	//! Should be a multiple of SHMEM_CHUNK_SIZE for better coalescing
	int _override_system_per_block;

	public:
	integrator(const config &cfg) : Base(cfg) {
		int spb = cfg.optional("system_per_block",0);
		_override_system_per_block =(spb % SHMEM_CHUNK_SIZE == 0) ? spb : 0;
	}

	/////////////////////// Logical parallelization parameters ////////////////////////////////////
	 
	const int& override_system_per_block()const{
		return _override_system_per_block;
	}

	/**
	 * Calculate number of worker threads needed for each system.
	 *
	 * This is set to number of greater value between number of pairs and number of
	 * coordinate components total. For example if nbod is 4 then there are 6 pairs of bodies, but
	 * there are 4x3=12 components per thread.
	 *
	 * This function is used by other helper functions to determine
	 * block size and to decompose the thread ID value inside the kernel.
	 *
	 */
	template<class T>
	static GENERIC int thread_per_system(T compile_time_param){
		const int nbod = T::n;
		const int body_comp = nbod * 3;
		const int pair_count = nbod * (nbod - 1) / 2;
		return (body_comp>pair_count) ? body_comp : pair_count;
	}


	/**
	 * Helper Function: Logical amount of shared memory needed for force calculations per thread
	 * This value is calculated from number of bodies, if this is used in a kernel, the actual
	 * amount of shared memory allocated might be different.
	 *
	 * \todo This assumes one particular implementation of Gravitation class.  Rewrite so that different gravitaiton classes (and different integrators) can request different ammounts of shared memory
	 *
	 */
	template<class T>
	static GENERIC int shmem_per_system(T compile_time_param){
		const int nbod = T::n;
		const int pair_count = nbod * (nbod - 1) / 2;
		return pair_count * 3  * 2 * sizeof(double);
	}


};

	
}
}
}

