#include "hermite_gpu_bpt.h"
#include "meta.hpp"
#include "storage.hpp"


namespace swarm {
namespace hermite_gpu_bpt {

	const int MAX_NBODIES = 10;
		
	__constant__ ensemble gpu_hermite_ens;

	inline __device__ static double inner_product(const double a[3],const double b[3]){
		return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
	}


	/**!
	 * helper function for Gravitation that operates on each component
	 * it gets scalar part of acceleration as input and calculates one component of
	 * acceleration and jerk at a time
	 *
	 */
	template<int nbod>
	__device__ static void Gravitation_component(int c,int i
			,double dx[3],double dv[3],double scalar,double rv
			,double (&acc)[3][nbod],double (&jerk)[3][nbod]){

	}

	template<int nbod>
		struct shared_data {
			double acc_length[nbod];
			double jerk_length[nbod];
//			double acc[3][(nbod*(nbod-1)/2]
//			double jerk[3][(nbod*(nbod-1)/2]
		} ;

	/*! 
	 * templatized function object to calculate acceleration and jerk
	 * It updates accleration and jerk for one body: bodid. this function
	 * object is body of a n*n loop. so it should get called for every pair
	 *
	 */
	template<int nbod>
	struct Gravitation {
		ensemble::systemref& sys;
		shared_data<nbod> &shared;
		double &acc,&jerk;
		const int i;
		const int c;
		__device__ Gravitation(const int bodid,const int c,ensemble::systemref& sys,shared_data<nbod> &shared, double &acc, double &jerk)
			:sys(sys),shared(shared),i(bodid),acc(acc),jerk(jerk),c(c){
				acc = 0;
				jerk = 0;
			}

		__device__ void calc_pair(int j,double &ac,double& je)const{
			if(i != j){

				if( c == 0 ) {
					double dx[3] =  { sys[j].p(0)- sys[i].p(0),sys[j].p(1)- sys[i].p(1), sys[j].p(2)- sys[i].p(2) };
					double dv[3] =  { sys[j].v(0)- sys[i].v(0),sys[j].v(1)- sys[i].v(1), sys[j].v(2)- sys[i].v(2) };

					// computing scalar part of the acceleration
					double r2 =  dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2] ;
					double rinv = rsqrt(r2)  / r2;

					shared.jerk_length[i] =  inner_product(dx,dv) * 3. / r2;
					shared.acc_length[i] =  sys[j].mass() * rsqrt(r2) / r2;
				}
				__syncthreads();

				double dxc = sys[j].p(c)- sys[i].p(c);
				double dvc = sys[j].v(c)- sys[i].v(c);
				ac = dxc* shared.acc_length[i];
				je = (dvc - dxc * shared.jerk_length[i] ) * shared.acc_length[i];

			} else {
				ac = 0;
				je = 0;
			}

		}

		__device__ void operator()(int j)const{
			double ac,je;
			calc_pair(j,ac,je);
			acc += ac;
			jerk += je;
		}
		__device__ void compute(){
			//for(int j = 0; j < nbod; j++) (*this)(j);
			Unroller<0,nbod>::step(*this);
		}
	};



	__device__ void copy(double src[3],double des[3]){ 
		des[0] = src[0], des[1] = src[1] , des[2] = src[2]; 
	}


	template<int nbod>
	__global__ void gpu_hermite_bpt_integrator_kernel(double destination_time, double time_step){
		
		ensemble &ens = gpu_hermite_ens;

		// this kernel will process specific body of specific system
		// getting essential pointers
		int sysid = ((blockIdx.z * gridDim.y + blockIdx.y) * gridDim.x + blockIdx.x) * blockDim.y + threadIdx.y;
		int sysid_in_block = threadIdx.y;
		int thr = threadIdx.x;
		// Body number
		int b = thr / 3 ;
		// Component number
		int c = thr % 3 ;

		bool p = true;

		if(sysid >= ens.nsys() || b >= nbod || c >= 3) { return; }

		ensemble::systemref sys ( ens[sysid] );


/*
		int ii = nbod - 1 - thr / (nbod/2);
		int jj = thr % (nbod/2);
		int b, j 
		if (jj < ii) {
			b = jj;
			j = ii;
		} else {
			b = nbod - 1 - jj - nbod%2;
			j = nbod - 1 - ii - nbod%2 + 1;
		}
		bool p = b < nbod;
*/


		// shared memory allocation
		extern __shared__ char shared_mem[];
		char * shared_memory_pointer = shared_mem + sysid_in_block * sizeof(shared_data<nbod>);

		shared_data<nbod> (&shared) = *( (shared_data<nbod>*) shared_memory_pointer );

//		double mass[nbod];
		// local memory allocation

		double t_start = sys.time(), t = t_start;
		double t_end = min(t_start + destination_time,sys.time_end());


		// Load data into shared memory (cooperative load)

		// local information per component per body
		double pos = sys[b].p(c), vel = sys[b].v(c), acc, jerk;

		

		// Calculate acceleration and jerk
		Gravitation<nbod> gi(b,c,sys,shared,acc,jerk);
		gi.compute(); //for(int j = 0; j < nbod; j++) gi(j);


		while(t < t_end){
			for(int k = 0; k < 2; k++)
			{
				double h = min(time_step, t_end - t);
				double pos_old = pos, vel_old = vel, acc_old = acc,jerk_old = jerk;

				// these two variable determine how each half of pos/vel/acc/jerk arrays
				// are going to be used to avoid unnecessary copying.
				// Predict 
				if( p ) {
					pos = pos_old +  h*(vel_old+(h*0.5)*(acc+(h/3.)*jerk));
					vel = vel_old +  h*(acc+(h*0.5)*jerk);
				}

			

				// Do evaluation and correction two times (PEC2)
				for(int l = 0; l < 2; l++)
				{

					// Gravitational force calculation

					// Write positions to shared memory
					if( p ) {
						sys[b].p(c) = pos, sys[b].v(c) = vel;
					}
					__syncthreads();

					// Calculate acceleration and jerk using shared memory
					Gravitation<nbod> gi(b,c,sys,shared,acc,jerk);
					gi.compute(); //for(int j = 0; j < nbod; j++) gi(j);


					// Correct
					if( p ) {
						pos = pos_old + (h*0.5) * ( (vel_old + vel) 
								+ (h*7.0/30.)*( (acc_old-acc) + (h/7.) * (jerk_old+jerk)));
						vel = vel_old + (h*0.5) * ( (acc_old+acc) + (h/6.) * (jerk_old-jerk));
					}

				}
				t += h;
			}

			if( p ) {
				sys[b].p(c) = pos, sys[b].v(c) = vel;
			}
		
			/*debug_hook();

			if(b == 0) 
				sys.increase_stepcount();
*/
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
			const int thread_per_system = ens.nbod() * 3 ;
			const int pair_count = ens.nbod() * (ens.nbod() - 1) / 2;
			const int shmem_per_system = ens.nbod() * 2 * sizeof(double);
			/* Upload the kernel parameters */ 
			if(ens.last_integrator() != this) 
			{ 
				int system_per_block = threadsPerBlock / thread_per_system;
				threadDim.x = thread_per_system;
				threadDim.y = system_per_block;
				
				// 2: new and old copies
				// 4: {pos,vel,acc,jerk}
				// 3: x,y,z
				this->shared_memory_size = system_per_block * shmem_per_system ;

				ens.set_last_integrator(this); 
				configure_grid(gridDim,  thread_per_system * system_per_block , ens.nsys()); 
				cudaMemcpyToSymbol(gpu_hermite_ens,	&ens, sizeof(gpu_hermite_ens) ); 
				if(dT == 0.) { return; } 
			} 
			// flush CPU/GPU output logs
			log::flush(log::memory | log::if_full);

			this->dT = dT;
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
