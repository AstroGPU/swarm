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

#include "swarm/common.hpp"
#include "swarm/gpu/bppt.hpp"
#include "monitors/stop_on_ejection.hpp"
#include "monitors/stop_on_any_large_distance_or_close_encounter.hpp"
#include "monitors/log_time_interval.hpp"
#include "monitors/combine.hpp"


namespace swarm {

template<int W = SHMEM_CHUNK_SIZE>
struct DoubleCoalescedStruct {
	typedef double scalar_t;
	double _value[W];
	GENERIC double& value(){ return _value[0]; }
};

namespace gpu {
namespace bppt {

/*! GPU implementation of PEC2 Hermite integrator w/ adaptive time step
 * \ingroup integrators
 *
 */
template< template<class L> class Monitor >
class hermite_adap: public integrator {
	typedef integrator base;
	typedef Monitor<gpulog::device_log> monitor_t;
	typedef typename monitor_t::params mon_params_t;
	private:
	double _time_step_factor, _min_time_step;
	mon_params_t _mon_params;

	public:
	hermite_adap(const config& cfg): base(cfg),_time_step_factor(0.001),_min_time_step(0.), _mon_params(cfg) {
		_time_step_factor =  cfg.require("time_step_factor", 0.0);
		_min_time_step =  cfg.require("min_time_step", 0.0);
	}

        // WARNING: Is this override working right?
	static GENERIC int shmem_per_system(int nbod) {
	   int mem_base = base::shmem_per_system(nbod);
	   int mem_derived = (2*3*nbod + nbod)*sizeof(double);
	   return std::max(mem_base,mem_derived);
	}

        // WARNING: Is this override working right?
	int  shmemSize(){
	   int mem_base = base::shmemSize();
	   int mem_derived = (system_per_block() * 2*3*_hens.nbod() + system_per_block()*_hens.nbod())*sizeof(double);
	   return std::max(mem_base,mem_derived);
	}

	virtual void launch_integrator() {
		launch_templatized_integrator(this);
	}


	template<class T>
	__device__ double calc_adaptive_time_step(T compile_time_param, SystemSharedData<T>& shared, const double acc, const double jerk)
		{
		// Body number
		int b = thread_body_idx(T::n);
		// Component number
		int c = thread_component_idx(T::n);
		bool body_component_grid = (b < T::n) && (c < 3);
		bool first_thread_in_system = (thread_in_system() == 0);

		// Shared Memory pointer computations here are wrong
		// The correct shared memory pointer calculation is done in
		// system_shared_data_pointer(T) in bppt.hpp
		// The system_shared_data_pointer can compute
		// extra shared memory for us as long as we define our data
		// structure using CoalescedDataStructure pattern

		// shared memory for computing time step
   //            double (&shared_acc_jerk_sq)[2][T::n][3] = * (double (*)[2][T::n][3]) system_shmem_1;
   //            double (&shared_time_step_factor)[T::n] = * (double (*)[T::n]) system_shmem_2;

		// Put accelerations and jerks for each body and component into shared memory
		if( body_component_grid ) {
		    shared.gravitation[b][c].acc() = acc*acc;
		    shared.gravitation[b][c].jerk() = jerk*jerk;
	        }
		__syncthreads();
		// calculate sum of squares of each component for each body
		// store ratio in shared memory
		if( (b < T::n) && (c==0) )
		    {
		    double acc_mag_sq = shared.gravitation[b][0].acc()+shared.gravitation[b][1].acc()+shared.gravitation[b][2].acc();
		    double jerk_mag_sq = shared.gravitation[b][0].jerk()+shared.gravitation[b][1].jerk()+shared.gravitation[b][2].jerk();
		    shared.time_step_factor[b].value() = jerk_mag_sq/acc_mag_sq;
		    }
		__syncthreads();
		if( first_thread_in_system ) 		    
		{
			double tf = shared.time_step_factor[0].value();
		    for(int bb=1;bb<T::n;++bb)
		      	tf += shared.time_step_factor[bb].value();
		    shared.time_step_factor[0].value() = rsqrt(tf)*_time_step_factor+_min_time_step;
		}
		__syncthreads();
		return shared.time_step_factor[0].value();
		}		

		template<class T>
		struct SystemSharedData {
			typename Gravitation<T::n>::shared_data gravitation;
			DoubleCoalescedStruct<> time_step_factor[T::n];
		};

	static GENERIC int shmem_per_system(int nbod){
		return bppt::shmem_per_system(nbod) + nbod*sizeof(double);
	}
	int  shmemSize(){
		const int nbod = _hens.nbod();
		// Round up number of systems in a block to the next multiple of SHMEM_CHUNK_SIZE
		int spb = ((system_per_block()+SHMEM_CHUNK_SIZE-1)/SHMEM_CHUNK_SIZE) * SHMEM_CHUNK_SIZE;
		return spb *  shmem_per_system(nbod);
	}


template< class T> 
static GPUAPI void * system_shared_data_pointer(T compile_time_param) {
	extern __shared__ char shared_mem[];
	int b = sysid_in_block() / SHMEM_CHUNK_SIZE ;
	int i = sysid_in_block() % SHMEM_CHUNK_SIZE ;
	int idx = i * sizeof(double) 
		+ b * SHMEM_CHUNK_SIZE 
		* shmem_per_system(T::n);
	return &shared_mem[idx];
}


	template<class T>
	__device__ void kernel(T compile_time_param){

		if(sysid()>=_dens.nsys()) return;
		// References to Ensemble and Shared Memory
		ensemble::SystemRef sys = _dens[sysid()];
		SystemSharedData* shared_data = (SystemSharedData*) system_shared_data_pointer(compile_time_param);
		Gravitation<T::n> calcForces(sys, shared_data->gravitation );

		// Local variables
		const int nbod = T::n;
		// Body number
		int b = thread_body_idx(nbod);
		// Component number
		int c = thread_component_idx(nbod);
		int ij = thread_in_system();
		bool body_component_grid = (b < nbod) && (c < 3);
		bool first_thread_in_system = (thread_in_system() == 0);


		// local variables
		monitor_t montest(_mon_params,sys,*_log) ;

		// local information per component per body
		double pos = 0.0, vel = 0.0 , acc0 = 0.0, jerk0 = 0.0;
		if( body_component_grid )
			pos = sys[b][c].pos() , vel = sys[b][c].vel();


		////////// INTEGRATION //////////////////////

		// Calculate acceleration and jerk
		calcForces(ij,b,c,pos,vel,acc0,jerk0);

		for(int iter = 0 ; (iter < _max_iterations) && sys.is_active() ; iter ++ ) {
			// Since h is the time step that is used for the step it makes more sense to
			// to calculate time step and assign it to h
			double h = calc_adaptive_time_step(compile_time_param, shared_data ,acc0,jerk0);

			if( sys.time() + h > _destination_time ) {
				h = _destination_time - sys.time();
			}

			
			// Initial Evaluation
			///calcForces(ij,b,c,pos,vel,acc0,jerk0);

			// Predict 
			pos = pos +  h*(vel+(h*0.5)*(acc0+(h/3.0)*jerk0));
			vel = vel +  h*(acc0+(h*0.5)*jerk0);

			double pre_pos = pos, pre_vel = vel;

			double acc1,jerk1;
			{
				// Evaluation
				calcForces(ij,b,c,pos,vel,acc1,jerk1);
				
				// Correct
				pos = pre_pos + (0.1-0.25) * (acc0 - acc1) * h * h - 1.0/60.0 * ( 7.0 * jerk0 + 2.0 * jerk1 ) * h * h * h;
				vel = pre_vel + ( -0.5 ) * (acc0 - acc1 ) * h -  1.0/12.0 * ( 5.0 * jerk0 + jerk1 ) * h * h;
				//	TODO: Need to test w/ new expressions below
				//				pos = pre_pos + ( (0.1-0.25) * (acc0 - acc1) - 1.0/60.0 * ( 7.0 * jerk0 + 2.0 * jerk1 ) * h) * h * h;
				// vel = pre_vel + (( -0.5 ) * (acc0 - acc1 ) -  1.0/12.0 * ( 5.0 * jerk0 + jerk1 ) * h )* h ;
			}
			{
				// Evaluation
				calcForces(ij,b,c,pos,vel,acc1,jerk1);
				
				// Correct
				pos = pre_pos + (0.1-0.25) * (acc0 - acc1) * h * h - 1.0/60.0 * ( 7.0 * jerk0 + 2.0 * jerk1 ) * h * h * h;
				vel = pre_vel + ( -0.5 ) * (acc0 - acc1 ) * h -  1.0/12.0 * ( 5.0 * jerk0 + jerk1 ) * h * h;
				//	TODO: Need to test w/ new expressions below
				// pos = pre_pos + ((0.1-0.25) * (acc0 - acc1) - 1.0/60.0 * ( 7.0 * jerk0 + 2.0 * jerk1 ) * h )* h * h ;
				// vel = pre_vel + (( -0.5 ) * (acc0 - acc1 ) -  1.0/12.0 * ( 5.0 * jerk0 + jerk1 ) * h ) * h ;
			}
			acc0 = acc1, jerk0 = jerk1;

			// Finalize the step
			if( body_component_grid )
				sys[b][c].pos() = pos , sys[b][c].vel() = vel;
			if( first_thread_in_system ) 
				sys.time() += h;

			if( first_thread_in_system && sys.is_active() )  {
				montest();
				if( sys.time() >= _destination_time ) 
					sys.set_inactive();
			}

			__syncthreads();


		}

	}


};


// WARNING: EBF: commented out to test new stopper
//integrator_plugin_initializer<hermite_adap< stop_on_ejection > >
//	hermite_adap_plugin("hermite_adap");

integrator_plugin_initializer<hermite_adap< monitors::stop_on_any_large_distance_or_close_encounter > >
	hermite_adap_plugin("hermite_adap");



}
}
}
