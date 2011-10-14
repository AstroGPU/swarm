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

	virtual void launch_integrator() {
		launch_templatized_integrator(this);
	}


	template<class T>
	__device__ double calc_adaptive_time_step(T compile_time_param, const double acc, const double jerk)
		{
		// Body number
		int b = thread_body_idx(T::n);
		// Component number
		int c = thread_component_idx(T::n);
		bool body_component_grid = (b < T::n) && (c < 3);
		bool first_thread_in_system = (thread_in_system() == 0);

		// shared memory for computing time step
		extern __shared__ char shared_mem[];
		// Array of squared acceleration and squared jerk for each body and each component
                char*  system_shmem_1 = (shared_mem + sysid_in_block() * sizeof(double)*2*T::n*3 );
                double (&shared_acc_jerk_sq)[2][T::n][3] = * (double (*)[2][T::n][3]) system_shmem_1;
		// Array of ratios for each body, sum to go in body 0
                char*  system_shmem_2 = (shared_mem + system_per_block_gpu() * sizeof(double)*2*T::n*3 + sysid_in_block() * sizeof(double)*T::n);
                double (&shared_time_step_factor)[T::n] = * (double (*)[T::n]) system_shmem_2;

		// Put accelerations and jerks for each body and component into shared memory
		if( body_component_grid ) {
		    shared_acc_jerk_sq[0][b][c] = acc*acc;
		    shared_acc_jerk_sq[1][b][c] = jerk*jerk;
		    }
		__syncthreads();
		// calculate sum of squares of each component for each body
		// store ratio in shared memory
		if( (b < T::n) && (c==0) )
		    {
		    double acc_mag_sq = shared_acc_jerk_sq[0][b][0]+shared_acc_jerk_sq[0][b][1]+shared_acc_jerk_sq[0][b][2];
		    double jerk_mag_sq = shared_acc_jerk_sq[1][b][0]+shared_acc_jerk_sq[1][b][1]+shared_acc_jerk_sq[1][b][2];
		    shared_time_step_factor[b] = jerk_mag_sq/acc_mag_sq;
		    }
		__syncthreads();
		if( first_thread_in_system ) 		    
		    {
		    for(int bb=1;bb<T::n;++bb)
		      shared_time_step_factor[0] += shared_time_step_factor[bb];
		    shared_time_step_factor[0] = rsqrt(shared_time_step_factor[0])*_time_step_factor+_min_time_step;
		    }
		__syncthreads();
		return shared_time_step_factor[0];
		}		


	template<class T>
	__device__ void kernel(T compile_time_param){

		if(sysid()>=_dens.nsys()) return;
		// References to Ensemble and Shared Memory
		ensemble::SystemRef sys = _dens[sysid()];
		typedef typename Gravitation<T::n>::shared_data grav_t;
//		Gravitation<T::n> calcForces(sys,*( (grav_t*) system_shared_data_pointer(compile_time_param) ) );
		Gravitation<T::n> calcForces(sys,sysid_in_block());

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
		double time_step = calc_adaptive_time_step(compile_time_param,acc0,jerk0);

		for(int iter = 0 ; (iter < _max_iterations) && sys.is_active() ; iter ++ ) {
			double h = time_step;

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
			time_step = calc_adaptive_time_step(compile_time_param,acc0,jerk0);

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
