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

/*! Paramaters for EulerPropagator
 * \ingroup propagator_parameters
 *
 */
struct EulerPropagatorParams {
	double time_step;
	EulerPropagatorParams(const config& cfg){
		time_step = cfg.require("time_step", 0.0);
	}
};

/*! GPU implementation of euler propagator
 * \ingroup propagators
 *
 */
template<class T>
struct EulerPropagator {
	typedef EulerPropagatorParams params;

	params _params;


	// Runtime variables
	ensemble::SystemRef& sys;
	Gravitation<T::n>& calcForces;
//	GravClass& calcForces;
	int b;
	int c;
	int ij;
	bool body_component_grid;
	bool first_thread_in_system;
	double max_timestep;


	GPUAPI EulerPropagator(const params& p,ensemble::SystemRef& s,
			Gravitation<T::n>& calc)
//			GravClass& calc)
		:_params(p),sys(s),calcForces(calc){}

	GPUAPI void init()  { }

	GPUAPI void shutdown() { }

        GPUAPI void convert_internal_to_std_coord() {} 
        GPUAPI void convert_std_to_internal_coord() {}

	GPUAPI void advance(){
		double h = min(_params.time_step, max_timestep);
		double pos = 0.0, vel = 0.0;
		double acc = 0.0, jerk = 0.0;
		const double third = 1.0/3.0;
		
		if( body_component_grid )
			pos = sys[b][c].pos() , vel = sys[b][c].vel();


		calcForces(ij,b,c,pos,vel,acc,jerk);
		// Integratore
		pos = pos +  h*(vel+(h*0.5)*(acc+(h*third)*jerk));
		vel = vel +  h*(acc+(h*0.5)*jerk);


		// Finalize the step
		if( body_component_grid )
			sys[b][c].pos() = pos , sys[b][c].vel() = vel;
		if( first_thread_in_system ) 
			sys.time() += h;
	}
};

typedef gpulog::device_log L;
using namespace monitors;

integrator_plugin_initializer< generic< EulerPropagator, stop_on_ejection<L> > >
	euler_prop_plugin("euler"
			,"This is the integrator based on euler propagator");

}
}
}
