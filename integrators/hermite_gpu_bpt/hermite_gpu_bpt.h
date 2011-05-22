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

/*! \file hermite_gpu_bpt.h
 *  \brief declares gpu_hermite_bpt_integrator
 *
 *  Note that while this clas derivers from integrator, it does not use gpu_generic_integrator
 */
#pragma once

#include <cuda_runtime_api.h>
#include "swarm.h"

namespace swarm {
namespace hermite_gpu_bpt {

/*!
 * \brief gpu_hermite_bpt_integrator class
 * computing only in double
 *
 */
class gpu_hermite_bpt_integrator : public integrator
{
	private:
	//// Variables
	ensemble* _gpu_ens;
	ensemble* _ens;

	//// Launch Variables
	int _threads_per_block;
	double _time_step;

	public:
	/*!
	 * \brief Constructor for hermite gpu integrator
	 *
	 * @param[in] cfg configuration class
	 */
	gpu_hermite_bpt_integrator(const config &cfg);

	~gpu_hermite_bpt_integrator() { if(_gpu_ens) cudaFree(_gpu_ens); }

	template<int nbod>
		void launch_template(const double& destination_time);

	/*!
	 * \brief host function to invoke a kernel (double precision) 
	 *
	 * @param[in,out] ens gpu_ensemble for data communication
	 * @param[in] dT destination time 
	 */
	void integrate(gpu_ensemble &ens, double dT){
		/* Upload ensemble */ 
		if(ens.last_integrator() != this) 
		{ 
			ens.set_last_integrator(this); 
			load_ensemble(ens);
		}

		launch_integrator(dT);
	}

	void load_ensemble(gpu_ensemble& ens){
		_ens = &ens;
		if(_gpu_ens)
			cudaFree(_gpu_ens);
		cudaMalloc(&_gpu_ens,sizeof(gpu_ensemble));
		cudaMemcpy(_gpu_ens, _ens, sizeof(gpu_ensemble),cudaMemcpyHostToDevice ); 
	}

	void launch_integrator(const double& destination_time);

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


} // end namespace hermite_gpu_bpt
} // end namespace swarm

