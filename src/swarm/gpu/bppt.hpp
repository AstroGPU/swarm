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
		* Impl::shmem_per_system(T::n);
	return &shared_mem[idx];
}

template<int W = SHMEM_CHUNK_SIZE>
struct DoubleCoalescedStruct {
	typedef double scalar_t;
	double _value[W];
	GENERIC double& value(){ return _value[0]; }
};


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

	///////////////////////// CUDA Grid Parameters ///////////////////////////////////////////

	/**
	 * Calculate the grid dimensions for CUDA kernel launch
	 * The generic kernel launcher calls this function to figure
	 * out the grid dimentions at the time of launch. This function
	 * can be overriden by descendants if a different grid is needed
	 */
	virtual dim3 gridDim(){
		const int nblocks = ( _hens.nsys() + system_per_block() - 1 ) / system_per_block();
		dim3 gD;
		gD.z = 1;
		find_best_factorization(gD.x,gD.y,nblocks);
		return gD;
	}

	/**
	 * Calculate the block dimensions for CUDA kernel launch
	 * The generic kernel launcher calls this function to figure
	 * out the grid dimentions at the time of launch. This function
	 * can be overriden by descendants if a different block size is needed
	 * 
	 * However, for most purposes, one may only need to change the 
	 * system_per_block and thread_per_system.
	 */
	virtual dim3 threadDim(){
		dim3 tD;
		tD.x = system_per_block();
		tD.y = thread_per_system();
		return tD;
	}

	/**
	 * Calculate amount of shared memory (used by generic kernel launcher)
	 * This method is used by generic kernel launcher at the time of launch to
	 * allocate dynamic shared memory. 
	 *
	 * An application that needs a different amount of shared memory can 
	 * override this function and return a different value
	 */
	int  shmemSize(){
		const int nbod = _hens.nbod();
		return system_per_block() * shmem_per_system(nbod);
	}


	/////////////////////// Logical parallelization parameters ////////////////////////////////////
	 
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
	virtual int thread_per_system() {
		const int nbod = _hens.nbod();
		const int body_comp = nbod * 3;
		const int pair_count = nbod * (nbod - 1) / 2;
		return std::max( body_comp, pair_count) ;
	}

	/**
	 * Calculate the number of systems that are fit into a block
	 *
	 * Although theoretically this can be set to 1, for performance reasons it is 
	 * better to set this value to something equal to the 
	 * number of memory banks in the device. This value is directly related to SHMEM_CHUNK_SIZE and
	 * ENSEMBLE_CHUNK_SIZE. For optimal performance all these constants should
	 * be equal to the number of memory banks in the device.
	 *
	 * Theoretical analysis: The load instruction is executed per warp. In Fermi
	 * architecture, a warp consists of 32 threads (CUDA 2.x) with consecutive IDs. A load
	 * instruction in a warp triggers 32 loads sent to the memory controller.
	 * Memory controller has a number of banks that can handle loads independently.
	 * According to CUDA C Programming Guide, there are no memory bank conflicts
	 * for accessing array of 64-bit valuse in devices with compute capabality of 2.x.
	 *
	 * It is essential to set system_per_block, SHMEM_CHUNK_SIZE and ENSEMBLE_CHUNK_SIZE to 
	 * the warp size for optimal performance. Lower values that are a power of 2 will result
	 * in some bank conflicts; the lower the value, higher is the chance of bank conflict.
	 */
	virtual int  system_per_block() {
		if(_override_system_per_block != 0)
			return _override_system_per_block;
		else {
			int tpc = thread_per_system();
			int cs = SHMEM_CHUNK_SIZE;
			int shm =  shmem_per_system(_hens.nbod());
			return optimized_system_per_block(cs, tpc, shm);
		}
	}

	/**
	 * Helper Function: Logical amount of shared memory needed for force calculations per thread
	 * This value is calculated from number of bodies, if this is used in a kernel, the actual
	 * amount of shared memory allocated might be different.
	 *
	 * \todo This assumes one particular implementation of Gravitation class.  Rewrite so that different gravitaiton classes (and different integrators) can request different ammounts of shared memory
	 *
	 */
	static GENERIC int shmem_per_system(int nbod) {
		const int pair_count = nbod * (nbod - 1) / 2;
		return pair_count * 3  * 2 * sizeof(double);
	}


};

	
}
}
}

