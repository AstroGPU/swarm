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

//// These functions are used inside kernels for consistent interpretation of thread and shared memory references

GPUAPI int sysid(){
	return ((blockIdx.z * gridDim.y + blockIdx.y) * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;
}
GPUAPI int sysid_in_block(){
	return threadIdx.x;
}
GPUAPI int thread_in_system() {
	return threadIdx.y;
}


GPUAPI int system_per_block_gpu() {
    return blockDim.x;
  };

GPUAPI int thread_component_idx(int nbod) {
	return thread_in_system() / nbod;
}
GPUAPI int thread_body_idx(int nbod) {
	return thread_in_system() % nbod;
}

GENERIC int shmem_per_system(int nbod) {
	const int pair_count = nbod * (nbod - 1) / 2;
	return pair_count * 3  * 2 * sizeof(double);
}

template< class T> 
GPUAPI void * system_shared_data_pointer(T compile_time_param) {
	extern __shared__ char shared_mem[];
	int b = sysid_in_block() / SHMEM_CHUCK_SIZE ;
	int i = sysid_in_block() % SHMEM_CHUCK_SIZE ;
	int idx = i * sizeof(double) 
		+ b * SHMEM_CHUCK_SIZE 
		* shmem_per_system(T::n);
	return &shared_mem[idx];
}

class integrator : public gpu::integrator  {
	typedef gpu::integrator Base;
	protected:
	//// Launch Variables
	int _system_per_block;

	public:
	integrator(const config &cfg) : Base(cfg) {
		// TODO: We may adjust systems_per_block to a lower value at launch time
		_system_per_block = cfg.optional("systems_per_block",16);
	}


	dim3 gridDim(){
		const int nblocks = ( _hens.nsys() + system_per_block() - 1 ) / system_per_block();
		dim3 gD;
		gD.z = 1;
		find_best_factorization(gD.x,gD.y,nblocks);
		return gD;
	}

	int thread_per_system() {
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

	int  shmemSize(){
		const int nbod = _hens.nbod();
		return system_per_block() * shmem_per_system(nbod);
	}


};

	
}
}
}

