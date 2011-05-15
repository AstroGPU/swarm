#include "hermite_gpu_bpt.h"
#include "static_accjerk.hpp"


namespace swarm {
namespace hermite_gpu_bpt {

	__constant__ ensemble gpu_hermite_ens;

	template<int nbod>
	__global__ void gpu_hermite_bpt_integrator_kernel(double destination_time, double time_step){
		
		ensemble &ens = gpu_hermite_ens;

		// this kernel will process specific body of specific system
		// getting essential pointers
		int sysid = ((blockIdx.z * gridDim.y + blockIdx.y) * gridDim.x + blockIdx.x) * blockDim.y + threadIdx.y;
		int sysid_in_block = threadIdx.y;
		int thr = threadIdx.x;

		if(sysid >= ens.nsys()) { return; }
		ensemble::systemref sys ( ens[sysid] );


		// Body/Component Grid
		// Body number
		int b = thr / 3 ;
		// Component number
		int c = thr % 3 ;
		bool body_component_grid = b < nbod;


		// i,j pairs Grid
		int ij = thr;


		// shared memory allocation
		extern __shared__ char shared_mem[];
		char * shared_memory_pointer = shared_mem + sysid_in_block * sizeof(Gravitation<nbod>::shared_data);

		struct Gravitation<nbod>::shared_data (&shared) = *( (struct Gravitation<nbod>::shared_data*) shared_memory_pointer );

//		double mass[nbod];
		// local memory allocation

		double t_start = sys.time(), t = t_start;
		double t_end = min(t_start + destination_time,sys.time_end());


		// local information per component per body
		double pos = 0, vel = 0 , acc = 0, jerk = 0;
		if( body_component_grid )
			pos = sys[b].p(c), vel = sys[b].v(c);

		

		// Calculate acceleration and jerk
		Gravitation<nbod> gi(sys,shared);
		gi(ij,b,c,acc,jerk);

		while(t < t_end){
			for(int k = 0; k < 2; k++)
			{
				double h = min(time_step, t_end - t);
				double pos_old = pos, vel_old = vel, acc_old = acc,jerk_old = jerk;

				// these two variable determine how each half of pos/vel/acc/jerk arrays
				// are going to be used to avoid unnecessary copying.
				// Predict 
				pos = pos_old +  h*(vel_old+(h*0.5)*(acc+(h/3.)*jerk));
				vel = vel_old +  h*(acc+(h*0.5)*jerk);

			

				// Do evaluation and correction two times (PEC2)
				for(int l = 0; l < 2; l++)
				{

					// Gravitational force calculation

					// Write positions to shared memory
					if( body_component_grid )
						sys[b].p(c) = pos, sys[b].v(c) = vel;
					__syncthreads();

					// Calculate acceleration and jerk using shared memory
					gi(ij,b,c,acc,jerk);


					// Correct
					pos = pos_old + (h*0.5) * ( (vel_old + vel) 
								+ (h*7.0/30.)*( (acc_old-acc) + (h/7.) * (jerk_old+jerk)));
					vel = vel_old + (h*0.5) * ( (acc_old+acc) + (h/6.) * (jerk_old-jerk));

				}
				t += h;
			}

			if( body_component_grid )
				sys[b].p(c) = pos, sys[b].v(c) = vel;
		
			if(log::needs_output(ens, t, sysid))
			{
				if(thr == 0) {
					sys.set_time(t);
					log::output_system(dlog, ens, t, sysid);
				}
			}

		}

		if(thr == 0) 
			sys.set_time(t);

	}


	template<int nbod>
		void gpu_hermite_bpt_integrator::integrate_internal(gpu_ensemble &ens, double dT)
		{
			if(ens.nbod() != nbod) return;

			dim3 gridDim;
			dim3 threadDim;

			const int body_comp = nbod * 3;
			const int pair_count = nbod * (nbod - 1) / 2;
			const int thread_per_system = max( body_comp, pair_count) ;
			const int shmem_per_system = pair_count * 3  * 2 * sizeof(double);
			const int system_per_block = threadsPerBlock / thread_per_system;
			const int shared_memory_size = system_per_block * shmem_per_system ;

			threadDim.x = thread_per_system;
			threadDim.y = system_per_block;
			configure_grid(gridDim,  thread_per_system * system_per_block , ens.nsys()); 

			// flush CPU/GPU output logs
			log::flush(log::memory | log::if_full);

			gpu_hermite_bpt_integrator_kernel<nbod><<<gridDim, threadDim, shared_memory_size>>>(dT, h);


			// flush CPU/GPU output logs
			log::flush(log::memory);

		}


	/*!
	 * \brief host function to invoke a kernel (double precision) 
	 *
	 * Currently maximum number of bodies is set to 10.
	 * In order to change, add if statement. 
	 * @param[in,out] ens gpu_ensemble for data communication
	 * @param[in] dT destination time 
	 */
		void gpu_hermite_bpt_integrator::integrate(gpu_ensemble &ens, double dT){
			/* Upload ensemble */ 
			if(ens.last_integrator() != this) 
			{ 
				ens.set_last_integrator(this); 
				cudaMemcpyToSymbol(gpu_hermite_ens,	&ens, sizeof(gpu_hermite_ens) ); 
				if(dT == 0.) { return; } 
			} 

			if(ens.nbod() <= 9){
				integrate_internal<3>(ens,dT);
				/*
				integrate_internal<4>(ens,dT);
				integrate_internal<5>(ens,dT);
				integrate_internal<6>(ens,dT);
				integrate_internal<7>(ens,dT);
				integrate_internal<8>(ens,dT);
				integrate_internal<9>(ens,dT);*/
			} else {
				// How do we get an error message out of here?
				ERROR("Invalid number of bodies. (only up to 10 bodies per system)");
				return;
			}
		}

}
}
