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

 // TODO: This probabaly belongs elsewhere
/*!
  \brief compute and store the number of systems that remain active

  NOTE: assumes not more than MAXTHREADSPERBLOCK threads per block
  NOTE: assumes a nthreads is a power of 2
  NOTE: assumes *nrunning = 0 on input
 @param[out] nrunning number of active systems
 @param[in] ens ensemble
*/
template< class implementation>
__global__ void count_running(int *nrunning, double Tstop, ensemble ens, implementation* integ )
{
        const int MAXTHREADSPERBLOCK = 256;

#if 0   // WARNING: This code doesn't appear to work
	// If you don't want to use shared memory for this
	// Direct counting (via atomicAdd) of the number of running threads. Slower
	// but requires no shared memory.
	// TODO: Test if this implementation is fast enough to be activated
	int sys = blockIdx.x * blockDim.x + threadIdx.x; // threadId();
	int running = (sys < ens.nsys()) && !(ens.flags(sys) & ensemble::INACTIVE || ens.time(sys) >= Tstop ) ? 1 : 0;
	if(running) { atomicAdd(nrunning, 1); /*printf("sys=%d running.\n", sys);*/ }
		__syncthreads();
		if(sys<3)	  printf("sys=%d nsys=%d flags=%d t=%g Tstop=%g running=%d nrunning=%d\n",sys,ens.nsys(),(ens.flags(sys) & ensemble::INACTIVE),ens.time(sys),Tstop,running,nrunning);
	return;
#else
	// We think it's ok to make this extern ?
	// TODO: replcae with dynamicly allocated shared memory 
	__shared__ int running[MAXTHREADSPERBLOCK];	// takes up 1k of shared memory (for MAXTHREADSPERBLOCK=256)

	int tpb = blockDim.x; //threadsPerBlock();
	int sys = blockIdx.x * blockDim.x + threadIdx.x; // threadId();
	const int widx = sys % tpb;
	running[widx] = (sys < ens.nsys()) && !(ens.flags(sys) & ensemble::INACTIVE || ens.time(sys) >= Tstop) ? 1 : 0;

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
#if 0
	__syncthreads();
	if(widx==0)
	  {
	  printf("sys=%d nsys=%d flags=%d t=%g Tstop=%g running=%d  ",sys,ens.nsys(),(ens.flags(sys) & ensemble::INACTIVE),ens.time(sys),Tstop,running[widx]);
	printf("nrunning=%d\n",nrunning[0]);	
	  }
#endif
#endif
}

class integrator : public swarm::integrator {

	protected:

	//// Launch Variables
	int _threads_per_block;
        unsigned int _max_itterations_per_kernel_call;
        unsigned int _max_kernel_calls_per_integrate_call;
	double _destination_time;

        /// integrator return type
        struct retval_t
	{
	  int nrunning;
	};

	public:
	integrator(const config &cfg) {
	  //		_threads_per_block = cfg.count("threads per block") ? atoi(cfg.at("threads per block").c_str()) : 128;
		_threads_per_block = cfg.count("blocksize") ? atoi(cfg.at("blocksize").c_str()) : 128;
		// TODO: Check that there are at least as many threads per block as threads per system;  For some reason adding assert causes segfaults
		//		assert(_threads_per_block>=thread_per_system());
		_max_itterations_per_kernel_call = cfg.count("max itterations per kernel call") ? atoi(cfg.at("max itterations per kernel call").c_str()) : 100000;
		_max_kernel_calls_per_integrate_call = cfg.count("max kernel calls per integrate call") ? atoi(cfg.at("max kernel calls per integrate call").c_str()) : 1000;

	}

	~integrator() { if(_gpu_ens) cudaFree(_gpu_ens); }

 // TODO: This probabaly belongs elsewhere
	int count_unfinished(gpu_ensemble &ens, double end_time )
           {
	     // initialize stop-time temporary array
	     //	     cuxDeviceAutoPtr<real_time> Tstop(ens.nsys(),end_time);
	     /// temp variable for return values (gpu pointer)
#if 1
	     cuxDeviceAutoPtr<int> nrunning_gpu(1);	
	     nrunning_gpu.memset(0);
	     int nsys = _ens->nsys();
	     int blocksize = 64;
	     int nblocks = (nsys/blocksize) + (nsys%blocksize!=0 ? 1 : 0);
	     dim3 gridDim(nblocks,1,1);
	     dim3 threadDim(blocksize,1,1);
	     count_running<<<gridDim, threadDim>>>(nrunning_gpu, end_time, ens, this);
	     int nrunning_cpu;
	     nrunning_gpu.download(&nrunning_cpu,1);
	     //	     std::cerr << "# end_time= " << end_time << " nrunning = " << nrunning_cpu << "\n";
	     return nrunning_cpu;
#else
	     return 1;
#endif
}

         void integrate(gpu_ensemble &ens, double stop_time ){
		/* Upload ensemble */ 
		if(ens.last_integrator() != this) 
		{ 
			ens.set_last_integrator(this); 
			load_ensemble(ens);

			// check if called integrate(ens,0.) to initialize integrator without integrating
			// TODO: Change the way dT and stop time are handeled
			//			if(dT == 0.) { return; }
			if(stop_time == 0.) { return; }
		}

		// TODO: Change the way dT and stop time are handeled
		_destination_time = stop_time;
		// advance_time_end_all(dT);
		// set_time_end_all(stop_time);

		// flush CPU/GPU output logs
		log::flush(log::memory | log::if_full);

		for(unsigned int l=0;l<_max_kernel_calls_per_integrate_call;++l)
		  {
		    launch_integrator();

		    // flush CPU/GPU output logs
		    log::flush(log::memory | log::if_full);

		    // TODO:  Need to replace this with something that checks how many systems are still trying to be integrated for dT
		    //		    unsigned int num_systems_active = 1;
		    unsigned int num_systems_active = count_unfinished(ens,stop_time);

		    //		    std::cerr << "# time = " << _ens->time(0) << "\n";
		    //		    std::cerr << "# l= " << l << " num_unfinished_sys = " << num_systems_active << "\n";
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



#if 0 // Would this be a good way to return values from device?  E.g., the number of unfinished systems, so we don't keep calling the kernel when there's nothing left to do
template<class implementation, class resultT>
void launch_templatized_integrator(implementation* integ, resultT* result){
#endif
template<class implementation>
void launch_templatized_integrator(implementation* integ){

	if(integ->get_ensemble()->nbod() <= 3){
		implementation* gpu_integ;
		cudaMalloc(&gpu_integ,sizeof(implementation));
		cudaMemcpy(gpu_integ,integ,sizeof(implementation),cudaMemcpyHostToDevice);

		params_t<3> gpu_integ_param;

#if 0  // Would this be a good way to return values from device?  E.g., the number of unfinished systems, so we don't keep calling the kernel when there's nothing left to do
		__device__ resultT d_result;
#endif
		launch_template(integ,gpu_integ,gpu_integ_param);
#if 0  // Would this be a good way to return values from device?  E.g., the number of unfinished systems, so we don't keep calling the kernel when there's nothing left to do
		if(sizeof(d_result)>0)
		  {
		    resultT h_result;
		    cudaMemcpyFromSymbol(&h_result,"d_result",sizeof(h_result),0,cudaMemcpyDeviceToHost);
		    *result = h_result;
		  }
#endif

		cudaFree(integ);
	} else {
		// How do we get an error message out of here?
		ERROR("Invalid number of bodies. (only up to 10 bodies per system)");
	}

}


	
}
}
