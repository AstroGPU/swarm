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

struct HermitePropagatorParams {
	double time_step;
	HermitePropagatorParams(const config& cfg){
		time_step = cfg.require("time step", 0.0);
	}
};

template<class T>
struct HermitePropagator {
	typedef HermitePropagatorParams params;

	params _params;


	// Runtime variables
	ensemble::SystemRef& sys;
	Gravitation<T::n>& calcForces;
	int b;
	int c;
	int ij;
	bool body_component_grid;
	bool first_thread_in_system;

	// Temporary variables
	double acc0, jerk0;

	GPUAPI HermitePropagator(const params& p,ensemble::SystemRef& s,
			Gravitation<T::n>& calc)
		:_params(p),sys(s),calcForces(calc){}

	GPUAPI void init() {
		double pos = 0, vel = 0;
		if( body_component_grid )
			pos = sys[b][c].pos() , vel = sys[b][c].vel();
		// Calculate acceleration and jerk
		calcForces(ij,b,c,pos,vel,acc0,jerk0);
	}

	GPUAPI void shutdown() {
	}

	GPUAPI void advance(){
		double h = _params.time_step;
		double pos = 0, vel = 0;
		if( body_component_grid )
			pos = sys[b][c].pos() , vel = sys[b][c].vel();

		// Predict 
		pos = pos +  h*(vel+(h*0.5)*(acc0+(h/3.)*jerk0));
		vel = vel +  h*(acc0+(h*0.5)*jerk0);

		double pre_pos = pos, pre_vel = vel;

		double acc1,jerk1;
		{
			// Evaluation
			calcForces(ij,b,c,pos,vel,acc1,jerk1);

			// Correct
			pos = pre_pos + (.1-.25) * (acc0 - acc1) * h * h - 1/60.0 * ( 7 * jerk0 + 2 * jerk1 ) * h * h * h;
			vel = pre_vel + ( -.5 ) * (acc0 - acc1 ) * h -  1/12.0 * ( 5 * jerk0 + jerk1 ) * h * h;
		}
		{
			// Evaluation
			calcForces(ij,b,c,pos,vel,acc1,jerk1);

			// Correct
			pos = pre_pos + (.1-.25) * (acc0 - acc1) * h * h - 1/60.0 * ( 7 * jerk0 + 2 * jerk1 ) * h * h * h;
			vel = pre_vel + ( -.5 ) * (acc0 - acc1 ) * h -  1/12.0 * ( 5 * jerk0 + jerk1 ) * h * h;
		}
		acc0 = acc1, jerk0 = jerk1;

		// Finalize the step
		if( body_component_grid )
			sys[b][c].pos() = pos , sys[b][c].vel() = vel;
		if( first_thread_in_system ) 
			sys.time() += h;
	}
};

integrator_plugin_initializer< generic< HermitePropagator, stop_on_ejection > >
	hermite_prop_plugin("hermite_prop"
			,"This is the integrator based on hermite propagator");

}
}
}
