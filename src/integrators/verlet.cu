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

#include "swarm/bppt.hpp"
#include "swarm/helpers.hpp"
#include "swarm/gravitation.hpp"
#include "monitors/stop_on_ejection.hpp"


namespace swarm {

namespace gpu {
namespace bppt {

template< template<class L> class Stopper >
class verlet: public integrator {
	typedef integrator base;
	typedef  typename Stopper<gpulog::device_log>::params stop_params_t;
	private:
	double _time_step;
	int _iteration_count;
	stop_params_t _stop_params;

	public:
	verlet(const config& cfg): base(cfg),_time_step(0.001), _stop_params(cfg) {
		if(!cfg.count("time step")) ERROR("Integrator gpu_verlet requires a timestep ('time step' keyword in the config file).");
		_time_step = atof(cfg.at("time step").c_str());
	}

	virtual void launch_integrator() {
		_iteration_count = _destination_time / _time_step;
		launch_templatized_integrator(this);
	}


	template<class T>
	__device__ void kernel(T a){

		if(sysid()>=_dens.nsys()) return;
		// References to Ensemble and Shared Memory
		ensemble::SystemRef sys = _dens[sysid()];
		typedef typename Gravitation<T::n>::shared_data grav_t;
		Gravitation<T::n> calcForces(sys,*( (grav_t*) system_shared_data_pointer(a) ) );

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
		Stopper<gpulog::device_log> stoptest(_stop_params,sys,*_log) ;


		// local information per component per body
		double pos = 0, vel = 0 , acc0 = 0, jerk0 = 0;
		if( body_component_grid )
			pos = sys[b][c].pos() , vel = sys[b][c].vel();

		double &t = sys.time();

		////////// INTEGRATION //////////////////////

		// Calculate acceleration and jerk
		calcForces(ij,b,c,pos,vel,acc0,jerk0);

		for(int iter = 0 ; (iter < _iteration_count) && sys.active() ; iter ++ ) {
			double h = _time_step;

			///////// INTEGRATION STEP /////////////////

			double h_first_half = h / 2 ;

			// First half step for positions
			pos = pos + h_first_half * vel;

			// Calculate acceleration in the middle
			double acc = calcForces.acc(ij,b,c,pos,vel);

			// First half step for velocities
			vel = vel + h_first_half * acc;

			// Step time to the middle of the step
			t  +=  h_first_half;

			// TODO: change half time step based on acc
			double h_second_half = h_first_half;

			// Second half step for positions and velocities
			vel = vel + h_second_half * acc;
			pos = pos + h_second_half * vel;

			// Step time to the end of the step
			t  +=  h_second_half;

			//////////////// END of Integration Step /////////////////


			// Finalize the step
			if( body_component_grid )
				sys[b][c].pos() = pos , sys[b][c].vel() = vel;
			if( first_thread_in_system ) 
				t += h;

			if( first_thread_in_system ) 
				sys.active() = ! stoptest() ;

			__syncthreads();


		}

	}


};

integrator_plugin_initializer<verlet< stop_on_ejection > >
	verlet_plugin("verlet");

}
}
}
