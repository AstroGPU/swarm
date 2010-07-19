#include "hermite_gpu_bpt.h"
#include "meta.hpp"

namespace swarm {
namespace hermite_gpu_bpt {

		
	__constant__ ensemble gpu_hermite_ens;

	inline __device__ static double inner_product(const double a[3],const double b[3]){
		return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
	}


	/**!
	 * helper function for accjerk_updater that operates on each component
	 * it gets scalar part of acceleration as input and calculates one component of
	 * acceleration and jerk at a time
	 *
	 */
	template<int nbod>
	__device__ static void accjerk_updater_component(int c,int i
			,double dx[3],double dv[3],double scalar,double rv
			,double (&acc)[3][nbod],double (&jerk)[3][nbod]){
		acc[c][i] += dx[c]* scalar;
		jerk[c][i] += (dv[c] - dx[c] * rv) * scalar;

	}

	/*! 
	 * templatized function object to calculate acceleration and jerk
	 * It updates accleration and jerk for one body: bodid. this function
	 * object is body of a n*n loop. so it should get called for every pair
	 *
	 */
	template<int nbod>
	struct accjerk_updater {
		ensemble::systemref& sysref;
		const double (&pos)[3][nbod],(&vel)[3][nbod];
		double (&acc)[3][nbod], (&jerk)[3][nbod];
		const int i;
		__device__ accjerk_updater(const int bodid,ensemble::systemref& sysref,const double (&pos)[3][nbod],const double (&vel)[3][nbod], double (&acc)[3][nbod], double (&jerk)[3][nbod])
			:sysref(sysref),pos(pos),vel(vel),acc(acc),jerk(jerk),i(bodid){
				acc[0][i] = acc[1][i] = acc[2][i] = 0.0;
				jerk[0][i] = jerk[1][i] = jerk[2][i] = 0.0;
			}
		__device__ void operator()(int j)const{
			if(i != j){

				double dx[3] =  { pos[0][j]-pos[0][i],pos[1][j]-pos[1][i],pos[2][j]-pos[2][i]};
				double dv[3] =  { vel[0][j]-vel[0][i],vel[1][j]-vel[1][i],vel[2][j]-vel[2][i]};

				// computing scalar part of the acceleration
				double r2 =  dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2] ;
				double rv =  inner_product(dx,dv) * 3 / r2;
				double rinv = rsqrt(r2)  / r2;

				// vectorized part
				const double scalar_i = +rinv*sysref[j].mass();
				for(int c = 0; c < 3 ; c++)
					accjerk_updater_component<nbod>(c,i,dx,dv,scalar_i,rv,acc,jerk);

			}
		}
		__device__ void compute(){
			for(int j = 0; j < nbod; j++) (*this)(j);
			//Unroller<0,nbod>::step(*this);
		}
	};


	template<int nbod>
		__device__ static void predictor(int i,int c,const int& s, const int& d
				,double (&pos)[2][3][nbod],double (&vel)[2][3][nbod],double (&acc)[2][3][nbod],double (&jerk)[2][3][nbod]
				,double h){
			pos[d][c][i] = pos[s][c][i] +  h*(vel[s][c][i]+(h/2)*(acc[s][c][i]+(h/3)*jerk[s][c][i]));
			vel[d][c][i] = vel[s][c][i] +  h*(acc[s][c][i]+(h/2)*jerk[s][c][i]);
		}

	template<int nbod>
	__device__ static void corrector(int i,const int& c,const int& s, const int& d
			,double (&pos)[2][3][nbod],double (&vel)[2][3][nbod],double (&acc)[2][3][nbod],double (&jerk)[2][3][nbod]
			,const double& h){
		pos[d][c][i] = pos[s][c][i] + (h/2) * ( (vel[s][c][i]+vel[d][c][i]) 
				+ (h*7.0/30)*( (acc[s][c][i]-acc[d][c][i]) + (h/7) * (jerk[s][c][i]+jerk[d][c][i])));
		vel[d][c][i] = vel[s][c][i] + (h/2) * ( (acc[s][c][i]+acc[d][c][i]) + (h/6) * (jerk[s][c][i]-jerk[d][c][i]));
	}

	__device__ void copy(double src[3],double des[3]){
		des[0] = src[0], des[1] = src[1] , des[2] = src[2];
	}

	template<int nbod>
	__device__ void load_to_shared(double pos[2][3][nbod],double vel[2][3][nbod],ensemble::systemref& sysref,int k,int c,int i){
		pos[k][c][i] = sysref[i].p(c), vel[k][c][i] = sysref[i].v(c);
	}

	template<int nbod>
	__device__ void store_from_shared(double pos[2][3][nbod],double vel[2][3][nbod],ensemble::systemref& sysref,int k,int c,int i){
		sysref[i].p(c) = pos[k][c][i] ,  sysref[i].v(c) = vel[k][c][i];
	}


	template<int nbod>
	__global__ void gpu_hermite_bpt_integrator_kernel(double destination_time, double time_step){
		
		// this kernel will process specific body of specific system
		// getting essential pointers
		int sysid = ((blockIdx.z * gridDim.y + blockIdx.y) * gridDim.x + blockIdx.x) * blockDim.y + threadIdx.y;
		int sysid_in_block = threadIdx.y;
		int bodid = threadIdx.x;
		ensemble &ens = gpu_hermite_ens;
		if(sysid >= ens.nsys() || bodid >= nbod) { return; }
		ensemble::systemref gsys ( ens[sysid] );


		// shared memory allocation
		extern __shared__ char shared_mem[];
		double (*shared_array)[2][3][nbod] = (double (*)[2][3][nbod]) shared_mem;
			
		// pointers to shared_memory
		double (&pos)[2][3][nbod] = shared_array[sysid_in_block*4], (&vel)[2][3][nbod] = shared_array[sysid_in_block*4+1];
		double (&acc)[2][3][nbod] = shared_array[sysid_in_block*4+2] , (&jerk)[2][3][nbod] = shared_array[sysid_in_block*4+3];

//		double mass[nbod];
		// local memory allocation

		double t_start = gsys.time(), t = t_start;
		double t_end = min(t_start + destination_time,gsys.time_end());

		// Load data into shared memory (cooperative load)
		load_to_shared<nbod>(pos,vel,gsys,0,0,bodid);
		load_to_shared<nbod>(pos,vel,gsys,0,1,bodid);
		load_to_shared<nbod>(pos,vel,gsys,0,2,bodid);
//		#pragma unroll
//		for(int i = 0; i < nbod; i++) mass[i] = gsys[i].mass();

		
		__syncthreads(); // load should complete before calculating acceleration and jerk

		// Calculate acceleration and jerk
		accjerk_updater<nbod> accjerk_updater_instance(bodid,gsys,pos[0],vel[0],acc[0],jerk[0]);
		accjerk_updater_instance.compute();

		while(t < t_end){
			for(int k = 0; k < 2; k++)
			{
				double h = min(time_step, t_end - t);
				// these two variable determine how each half of pos/vel/acc/jerk arrays
				// are going to be used to avoid unnecessary copying.
				const int s = k, d = 1-k; 
				// Predict 
				for(int c = 0; c < 3; c++)
					predictor<nbod>(bodid,c,s,d,pos,vel,acc,jerk,h);


				// Do evaluation and correction two times (PEC2)
				for(int l = 0; l < 2; l++)
				{
					__syncthreads();
					// Calculate acceleration and jerk
					accjerk_updater<nbod> accjerk_updater_instance(bodid,gsys,pos[d],vel[d],acc[d],jerk[d]);
					accjerk_updater_instance.compute();
					//Unroller<0,nbod>::step(accjerk_updater_instance);

					//__syncthreads(); // to prevent WAR. corrector updates pos/vel that accjerk_updater would read

					// Correct
					for(int c = 0; c < 3; c++)
						corrector<nbod>(bodid,c,s,d,pos,vel,acc,jerk,h);

				}
				t += h;
			}

			/*debug_hook();

			if(bodid == 0) 
				gsys.increase_stepcount();
*/
			if(log::needs_output(ens, t, sysid))
			{
				// Save pos/vel to global memory
				store_from_shared<nbod>(pos,vel,gsys,0,0,bodid);
				store_from_shared<nbod>(pos,vel,gsys,0,1,bodid);
				store_from_shared<nbod>(pos,vel,gsys,0,2,bodid);
				if(bodid == 0) {
					gsys.set_time(t);
					log::output_system(dlog, ens, t, sysid);
				}
			}

		}

		if(bodid == 0) 
			gsys.set_time(t);
		// Save pos/vel to global memory
		store_from_shared<nbod>(pos,vel,gsys,0,0,bodid);
		store_from_shared<nbod>(pos,vel,gsys,0,1,bodid);
		store_from_shared<nbod>(pos,vel,gsys,0,2,bodid);

	}

	//! Simple template Function Object to execute appropriate kernel at runtime
	template<int nbod>
		struct kernel_launcher {
			template<class P>
				static void choose(P p){
					double dT = p.dT;
					double h = p.h;
					gpu_hermite_bpt_integrator_kernel<nbod><<<p.gridDim, p.threadDim, p.shared_memory_size>>>(dT, h);
				}
		};



	inline int min_power_2(const int& x){
		int y;
		for(y = 1; y < x; y*=2);
		return y;
	}
	/*!
	 * \brief host function to invoke a kernel (double precision) 
	 *
	 * Currently maximum number of bodies is set to 10.
	 * In order to change, add if statement. 
	 * @param[in,out] ens gpu_ensemble for data communication
	 * @param[in] dT destination time 
	 */
		void gpu_hermite_bpt_integrator::integrate(gpu_ensemble &ens, double dT)
		{
			/* Upload the kernel parameters */ 
			if(ens.last_integrator() != this) 
			{ 
				int system_per_block = threadsPerBlock / ens.nbod();
				threadDim.x = ens.nbod();
				threadDim.y = system_per_block;
				
				// 2: new and old copies
				// 4: {pos,vel,acc,jerk}
				// 3: x,y,z
				this->shared_memory_size = system_per_block * ens.nbod() * 2 * 4 * 3 * sizeof(double);

				ens.set_last_integrator(this); 
				configure_grid(gridDim,  system_per_block , ens.nsys()); 
				cudaMemcpyToSymbol(gpu_hermite_ens,	&ens, sizeof(gpu_hermite_ens) ); 
				if(dT == 0.) { return; } 
			} 
			// flush CPU/GPU output logs
			log::flush(log::memory | log::if_full);

			this->dT = dT;
			const int MAX_NBODIES = 10;
			if(ens.nbod() <= MAX_NBODIES){
				choose<kernel_launcher,3,MAX_NBODIES,void>(ens.nbod(),*this);
			} else {
				// How do we get an error message out of here?
				ERROR("Invalid number of bodies. (only up to 10 bodies per system)");
				return;
			}
			// flush CPU/GPU output logs
			log::flush(log::memory);

		}

}
}
