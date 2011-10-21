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

/*! GPU implementation of PEC2 Hermite integrator
 * \ingroup integrators
 *
 */
template< template<class L> class Monitor >
class hermite: public integrator {
	typedef integrator base;
	typedef Monitor<gpulog::device_log> monitor_t;
	typedef typename monitor_t::params mon_params_t;
	private:
	double _time_step;
	mon_params_t _mon_params;

	public:
	hermite(const config& cfg): base(cfg),_time_step(0.001), _mon_params(cfg) {
		_time_step =  cfg.require("time_step", 0.0);
	}

	virtual void launch_integrator() {
		launch_templatized_integrator(this);
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
		bool first_thread_in_system = thread_in_system() == 0;


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
			double h = _time_step;

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

			if( first_thread_in_system  )  {
			    montest();
			    if( sys.time() >= _destination_time ) 
				sys.set_inactive();
			}

			__syncthreads();


		}

	}


};


// WARNING: EBF: commented out to test new stopper
//integrator_plugin_initializer<hermite< stop_on_ejection > >
//	hermite_plugin("hermite");

integrator_plugin_initializer<hermite< monitors::stop_on_any_large_distance_or_close_encounter > >
	hermite_plugin("hermite");

integrator_plugin_initializer<hermite< monitors::log_time_interval > >
	hermite_log_plugin("hermite_log");


}
}
}
