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
#include "datatypes.hpp"
#include "ensemble.hpp"


namespace swarm {
namespace hp {

class integrator  {

	protected:
	//// Launch Variables
	int _system_per_block;
	double _destination_time;
	hostEnsemble _hens;
	//TODO: use cux auto ptr to make sure we don't have memory leaks
	deviceEnsemble _dens;

	public:
	integrator(const config &cfg) {
		_system_per_block = cfg.count("blocksize") ? atoi(cfg.at("blocksize").c_str()) : 32;
	}

	ensemble& get_ensemble() {
		return _hens;
	}

	void load_ensemble(hostEnsemble& ens){
		_hens = ens;
		_dens = _hens.cloneTo<deviceEnsemble>();
	}
	
	void download_ensemble() {
		_dens.copyTo( _hens );
	}

	void set_duration(const double& duration) {
		_destination_time = duration;
	}

	void load_ensemble(cpu_ensemble& ens){
		_hens = hostEnsemble::create(ens.nbod(),ens.nsys());

		// Copy all
		for(int i= 0; i < ens.nsys() ; i++){
			for(int j = 0; j < ens.nbod(); j++) {
				_hens[i][j][0].pos() = ens[i][j].p(0);
				_hens[i][j][1].pos() = ens[i][j].p(1);
				_hens[i][j][2].pos() = ens[i][j].p(2);
				_hens[i][j][0].vel() = ens[i][j].v(0);
				_hens[i][j][1].vel() = ens[i][j].v(1);
				_hens[i][j][2].vel() = ens[i][j].v(2);
				_hens[i][j].mass() = ens[i][j].mass();
				_hens[i][j].time() = ens[i].time();
			}
		}

		load_ensemble(_hens);
	}


	void integrate(double dT){
		_destination_time = dT;

		// flush CPU/GPU output logs
		log::flush(log::memory | log::if_full);

		launch_integrator();

		// flush CPU/GPU output logs
		log::flush(log::memory);
	}


	virtual void launch_integrator() = 0;

	dim3 gridDim(){
		const int nblocks = ( _hens.nsys() + system_per_block() ) / system_per_block();
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

	static GPUAPI int shmem_per_system(int nbod) {
		const int pair_count = nbod * (nbod - 1) / 2;
		return pair_count * 3  * 2 * sizeof(double);
	}

};

inline __device__ int sysid(){
	return ((blockIdx.z * gridDim.y + blockIdx.y) * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;
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
	
}
}
