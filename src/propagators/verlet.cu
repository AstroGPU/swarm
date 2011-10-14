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

#include "swarm/swarmplugin.h"



namespace swarm {

namespace gpu {
namespace bppt {

struct VerletPropagatorParams {
	double time_step;
	VerletPropagatorParams(const config& cfg){
		time_step = cfg.require("time_step", 0.0);
	}
};

template<class T>
//template<class T, class GravClass>
struct VerletPropagator {
	typedef VerletPropagatorParams params;

	params _params;


	// Runtime variables
	ensemble::SystemRef& sys;
	Gravitation<T::n>& calcForces;
//	GravitationAccOnly<T::n>& calcForces;
//	GravClass& calcForces;
	int b;
	int c;
	int ij;
	bool body_component_grid;
	bool first_thread_in_system;
	double max_timestep;


	GPUAPI VerletPropagator(const params& p,ensemble::SystemRef& s,
			Gravitation<T::n>& calc)
//			GravitationAccOnly<T::n>& calc)
//			GravClass& calc)
		:_params(p),sys(s),calcForces(calc){}

	GPUAPI void init()  { }

	GPUAPI void shutdown() { }

	GPUAPI void advance(){
		double h = _params.time_step;
		double pos = 0.0, vel = 0.0;

		if( body_component_grid )
			pos = sys[b][c].pos() , vel = sys[b][c].vel();


			///////// INTEGRATION STEP /////////////////

			double h_first_half = 0.5 * h;

			// First half step for positions
			pos = pos + h_first_half * vel;

			// Calculate acceleration in the middle
			double acc = calcForces.acc(ij,b,c,pos,vel);

			// First half step for velocities
			vel = vel + h_first_half * acc;

			// TODO: change half time step based on acc
			double h_second_half = h_first_half;

			// Second half step for positions and velocities
			vel = vel + h_second_half * acc;
			pos = pos + h_second_half * vel;

			//////////////// END of Integration Step /////////////////


		// Finalize the step
		if( body_component_grid )
			sys[b][c].pos() = pos , sys[b][c].vel() = vel;
		if( first_thread_in_system ) 
			sys.time() += h_first_half + h_second_half;
	}
};

integrator_plugin_initializer< generic< VerletPropagator, monitors::stop_on_ejection > >
	verlet_prop_plugin("verlet"
			,"This is the integrator based on verlet propagator");

}
}
}
