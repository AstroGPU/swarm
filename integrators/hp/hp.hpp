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

#include <cuda_runtime_api.h>
#include "swarm.h"
#include "swarmlog.h"
#include "gravitation.hpp"

namespace swarm {
namespace hp {

class integrator : public swarm::integrator {

	protected:

	//// Launch Variables
	int _threads_per_block;
        unsigned int _max_itterations_per_kernel_call;
        unsigned int _max_kernel_calls_per_integrate_call;
	double _destination_time;

	public:
	integrator(const config &cfg) {
	  //		_threads_per_block = cfg.count("threads per block") ? atoi(cfg.at("threads per block").c_str()) : 128;
		_threads_per_block = cfg.count("blocksize") ? atoi(cfg.at("blocksize").c_str()) : 128;
		// TODO: Check that there are at least as many threads per block as threads per system;  For some reason adding assert causes segfaults
		//		assert(_threads_per_block>=thread_per_system());
		_max_itterations_per_kernel_call = cfg.count("max itterations per kernel call") ? atoi(cfg.at("max itterations per kernel call").c_str()) : 10000;
		// TODO: Once are able to check for how many systems are still active, increase the default value
		_max_kernel_calls_per_integrate_call = cfg.count("max kernel calls per integrate call") ? atoi(cfg.at("max kernel calls per integrate call").c_str()) : 10;

	}

	~integrator() { if(_gpu_ens) cudaFree(_gpu_ens); }

         void integrate(gpu_ensemble &ens, double dT ){
		/* Upload ensemble */ 
		if(ens.last_integrator() != this) 
		{ 
			ens.set_last_integrator(this); 
			load_ensemble(ens);

			// check if called integrate(ens,0.) to initialize integrator without integrating
			if(dT == 0.) { return; }
		}

		_destination_time = dT;

		// flush CPU/GPU output logs
		log::flush(log::memory | log::if_full);

		for(unsigned int l=0;l<_max_kernel_calls_per_integrate_call;++l)
		  {
		    launch_integrator();

		    // flush CPU/GPU output logs
		    log::flush(log::memory | log::if_full);

		    // TODO:  Need to replace this with something that checks how many systems are still trying to be integrated for dT
		    unsigned int num_systems_active = 1;
		    if(num_systems_active<=0)
		      break;
		  }

		// flush CPU/GPU output logs
		log::flush(log::memory);
	}


	virtual void launch_integrator() = 0;

	int thread_per_system() const {
		const int nbod = _ens->nbod();
		const int body_comp = nbod * 3;
		const int pair_count = nbod * (nbod - 1) / 2;
		const int threads_per_system = std::max( body_comp, pair_count) ;
		return threads_per_system;
	}

	int  system_per_block() const {
		return  _threads_per_block / thread_per_system();
	}

	dim3 gridDim()  const {
		const int threads_per_system = thread_per_system();
		const int systems_per_block = system_per_block();
		// TODO: Saleh, please checkk that my replacement is correct
		// const int nblocks = ( _ens->nsys() + systems_per_block ) / systems_per_block;
		const int nblocks = ( _ens->nsys() / systems_per_block ) +  ( (_ens->nsys() % systems_per_block != 0 ) ? 1 : 0 );



		dim3 gD;
		gD.z = 1;
		find_best_factorization(gD.x,gD.y,nblocks);
		return gD;
	}
	dim3 threadDim() const {
	        const int threads_per_system = thread_per_system();
	  	const int systems_per_block = system_per_block();

		dim3 tD;
		tD.x = threads_per_system;
		tD.y = systems_per_block;
		return tD;
	}

        int  shmemSize() const {
		const int nbod = _ens->nbod();
		const int systems_per_block = system_per_block();
		return systems_per_block * shmem_per_system(nbod);
	}

	static __device__ __host__ inline int shmem_per_system(int nbod) {
		const int pair_count = nbod * (nbod - 1) / 2;
		return pair_count * 3  * 2 * sizeof(double);
	}

};

template<int i>
struct params_t {
	const static int n = i;
};

template<class implementation,class T>
__global__ void generic_kernel(implementation* integ,T a) {
	integ->kernel(a);
}


template< class implementation, class T>
void launch_template(implementation* integ, implementation* gpu_integ, T a)
{
	if(integ->get_ensemble()->nbod() == T::n) 
		generic_kernel<<<integ->gridDim(), integ->threadDim(), integ->shmemSize() >>>(gpu_integ,a);

}

template<class implementation>
void launch_templatized_integrator(implementation* integ){

	if(integ->get_ensemble()->nbod() <= 3){
		implementation* gpu_integ;
		cudaMalloc(&gpu_integ,sizeof(implementation));
		cudaMemcpy(gpu_integ,integ,sizeof(implementation),cudaMemcpyHostToDevice);

		launch_template(integ,gpu_integ,params_t<3>());

		cudaFree(integ);
	} else {
		// How do we get an error message out of here?
		ERROR("Invalid number of bodies. (only up to 10 bodies per system)");
	}

}


	
}
}
