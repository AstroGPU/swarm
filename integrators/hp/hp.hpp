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

namespace swarm {
namespace hp {

class integrator : public swarm::integrator {

	protected:
	//// Launch Variables
	int _threads_per_block;

	public:
	integrator(const config &cfg) {
		_threads_per_block = cfg.count("threads per block") ? atoi(cfg.at("threads per block").c_str()) : 128;
	}

	~integrator() { if(_gpu_ens) cudaFree(_gpu_ens); }

	void integrate(gpu_ensemble &ens, double dT){
		/* Upload ensemble */ 
		if(ens.last_integrator() != this) 
		{ 
			ens.set_last_integrator(this); 
			load_ensemble(ens);
		}

		launch_integrator(dT);
	}


	virtual void launch_integrator(const double& destination_time) = 0;

	dim3 gridDim(){
		const int nbod = _ens->nbod();
		const int body_comp = nbod * 3;
		const int pair_count = nbod * (nbod - 1) / 2;
		const int thread_per_system = std::max( body_comp, pair_count) ;
		const int system_per_block = _threads_per_block / thread_per_system;
		const int nblocks = ( _ens->nsys() + system_per_block ) / system_per_block;

		dim3 gD;
		gD.z = 1;
		find_best_factorization(gD.x,gD.y,nblocks);
		return gD;
	}
	dim3 threadDim(){
		const int nbod = _ens->nbod();
		const int body_comp = nbod * 3;
		const int pair_count = nbod * (nbod - 1) / 2;
		const int thread_per_system = std::max( body_comp, pair_count) ;
		const int system_per_block = _threads_per_block / thread_per_system;

		dim3 tD;
		tD.x = thread_per_system;
		tD.y = system_per_block;
		return tD;
	}
	int  shmemSize(){
		const int nbod = _ens->nbod();
		const int body_comp = nbod * 3;
		const int pair_count = nbod * (nbod - 1) / 2;
		const int thread_per_system = std::max( body_comp, pair_count) ;
		const int shmem_per_system = pair_count * 3  * 2 * sizeof(double);
		const int system_per_block = _threads_per_block / thread_per_system;
		const int shared_memory_size = system_per_block * shmem_per_system ;
		return shared_memory_size;
	}

};
	
}
}
