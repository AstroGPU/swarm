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

			struct MidpointPropagatorParams {
				double time_step;
				MidpointPropagatorParams(const config& cfg){
					time_step = cfg.require("time step", 0.0);
				}
			};

			template<class T>
				struct MidpointPropagator {
					typedef MidpointPropagatorParams params;

					params _params;


					// Runtime variables
					ensemble::SystemRef& sys;
					Gravitation<T::n>& calcForces;
					int b;
					int c;
					int ij;
					bool body_component_grid;
					bool first_thread_in_system;


					GPUAPI MidpointPropagator(const params& p,ensemble::SystemRef& s,
							Gravitation<T::n>& calc)
						:_params(p),sys(s),calcForces(calc){}

					GPUAPI void init()  { }

					GPUAPI void shutdown() { }

					GPUAPI void advance(){
						double H = _params.time_step;
						double pos = 0, vel = 0;

						if( body_component_grid )
							pos = sys[b][c].pos() , vel = sys[b][c].vel();


						////////// INTEGRATION //////////////////////

						/// Modified midpoint method integrator with n substeps
						const int n = 4;
						double h = H / n;

						double p_i , p_im1, p_im2;
						double v_i,  v_im1, v_im2;
						double a_im1;

						// Step 0
						p_i = pos;
						v_i = vel;

						// Step 1
						p_im1 = p_i;
						v_im1 = v_i;

						a_im1 = calcForces.acc(ij,b,c,p_im1,v_im1);

						p_i = p_im1 + h * v_im1;
						v_i = v_im1 + h * a_im1;

						// Step 2 .. n
						for(int i = 2; i <= n; i++){
							p_im2 = p_im1;
							p_im1 = p_i;
							v_im2 = v_im1;
							v_im1 = v_i;

							a_im1 = calcForces.acc(ij,b,c,p_im1,v_im1);

							p_i = p_im2 + 2 * h * v_im1;
							v_i = v_im2 + 2 * h * a_im1;
						}
						double a_i = calcForces.acc(ij,b,c,p_i,v_i);

						pos = ( p_i + p_im1 + h * v_i ) / 2;
						vel = ( v_i + v_im1 + h * a_i ) / 2;


						// Finalize the step
						if( body_component_grid )
							sys[b][c].pos() = pos , sys[b][c].vel() = vel;
						if( first_thread_in_system ) 
							sys.time() += H;
					}
				};

			integrator_plugin_initializer< generic< MidpointPropagator, stop_on_ejection > >
				midpoint_prop_plugin("midpoint"
						,"This is the integrator based on midpoint propagator");

		}
	}
}
