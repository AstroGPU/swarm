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
#include "gravitation.hpp"


namespace swarm {
namespace gpu {

/*! Class of GPU integrators with a thread for each body-pair
 * \addtogroup integrators
 *
 */
namespace bppt {


inline __device__ int sysid(){
	return ((blockIdx.z * gridDim.y + blockIdx.y) * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;
	// TODO: Test whether we could use just
	//	return ((blockIdx.y) * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;
}
inline __device__ int sysid_in_block(){
	return threadIdx.x;
}
inline __device__ int thread_in_system() {
	return threadIdx.y;
}

inline __device__ int thread_component_idx(int nbod) {
	return thread_in_system() / nbod;
}
inline __device__ int thread_body_idx(int nbod) {
	return thread_in_system() % nbod;
}

class integrator : public gpu::integrator  {
	typedef gpu::integrator Base;
	protected:
	//// Launch Variables
        int _system_per_block;
        int _thread_per_system;

	public:
	integrator(const config &cfg) : Base(cfg) {
	// TODO: Should we allow integrator to provide its own default block size?
	//       Or is this about memory mangagment, so not integrator specific?
	// WARNING: "blocksize" or "block size"?
	_system_per_block = cfg.count("blocksize") ? atoi(cfg.at("blocksize").c_str()) : 32;
	// WARNING: Should the above be changed to match what's in utils.cpp:validate_configuration? (see below)
	//  int bs = cfg.optional("block size",SHMEM_WARPSIZE);

	_thread_per_system = calc_thread_per_system();
	}

	dim3 gridDim(){
		// TODO:  Is this one too many blocks if nsys%system_per_block==0?
		const int nblocks = ( _hens.nsys() + system_per_block() ) / system_per_block();
		dim3 gD;
		gD.z = 1;
		find_best_factorization(gD.x,gD.y,nblocks);
		return gD;
	}

	int calc_thread_per_system() {
		const int nbod = _hens.nbod();
		const int body_comp = nbod * 3;
		const int pair_count = nbod * (nbod - 1) / 2;
		return std::max( body_comp, pair_count) ;
	}

	dim3 threadDim(){
		dim3 tD;
		tD.x = system_per_block();
		tD.y = thread_per_system();
		return tD;
	}

	int  system_per_block() {
		return  _system_per_block ;
	}

	int  thread_per_system() {
		return  _thread_per_system ;
	}

  // TODO: Do we really need to provide a default version (that we override if we use a non-default Gravitation class)?
  // Or should we make integrator take a GravitationClass as a template parameter, so we will already know how much memory it will need?
	int  shmemSize(){
	  //  WARNING: Is this right?  
	  //   What if the number of systems per block is not a multiple of SHMEM_WARPSIZE?  
	  //   If we enforce that it is, then there's no problem
		const int nbod = _hens.nbod();
		return system_per_block() * shmem_per_system(nbod);
	}


  // Shouldn't this be moved into gravitation class?
	static GENERIC int shmem_per_system(int nbod) {
		const int pair_count = nbod * (nbod - 1) / 2;
		return pair_count * 3  * 2 * sizeof(double);
	}

  // Shouldn't this be moved into gravitation class?
	template< class T> 
	static __device__ void * system_shared_data_pointer(T compile_time_param) {
		extern __shared__ char shared_mem[];
		int b = sysid_in_block() / SHMEM_WARPSIZE ;
		int i = sysid_in_block() % SHMEM_WARPSIZE ;
		int idx = i * sizeof(double) 
			+ b * SHMEM_WARPSIZE 
			* shmem_per_system(T::n);
		return &shared_mem[idx];
	}

  // Shouldn't this be moved into gravitation class?
	template< class T> 
	__device__ void * unused_shared_data_pointer(T compile_time_param) {
		extern __shared__ char shared_mem[];
		int idx = system_per_block() * shmem_per_system(T::n);
		return &shared_mem[idx];
	}


};

	
}
}
}
