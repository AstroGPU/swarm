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

/*! \file midpoint.hpp
 *   \brief Defines \ref swarm::gpu::bppt::MidpointPropagator - the GPU 
 *          implementation of modified midpoint method propagator. 
 *
 */

#include "swarm/swarmplugin.h"

namespace swarm { namespace gpu { namespace bppt {

/*! Paramaters for MidpointPropagator
 * \ingroup propagator_parameters
 *
 */
struct MidpointPropagatorParams {
	double time_step;
	MidpointPropagatorParams(const config& cfg){
		time_step = cfg.require("time_step", 0.0);
	}
};

/*! GPU implementation of modified midpoint method propagator
 * \ingroup propagators
 *
 */
template<class T,class Gravitation>
struct MidpointPropagator {
	typedef MidpointPropagatorParams params;
	const static int nbod = T::n;

	params _params;


	// Runtime variables
	ensemble::SystemRef& sys;
	Gravitation& calcForces;
	int b;
	int c;
	int ij;
//	bool body_component_grid;
//	bool first_thread_in_system;
	double max_timestep;


	GPUAPI MidpointPropagator(const params& p,ensemble::SystemRef& s,
			Gravitation& calc)
		:_params(p),sys(s),calcForces(calc){}

	static GENERIC int thread_per_system(){
		return nbod * 3;
	}

	static GENERIC int shmem_per_system() {
		 return 0;
	}

	GPUAPI void init()  { }

	GPUAPI void shutdown() { }

        GPUAPI void convert_internal_to_std_coord() {} 
        GPUAPI void convert_std_to_internal_coord() {}

	__device__ bool is_in_body_component_grid()
//        { return body_component_grid; }	
        { return  ((b < T::n) && (c < 3)); }	

	__device__ bool is_in_body_component_grid_no_star()
//        { return ( body_component_grid && (b!=0) ); }	
        { return ( (b!=0) && (b < T::n) && (c < 3) ); }	

	__device__ bool is_first_thread_in_system()
//        { return first_thread_in_system; }	
        { return (thread_in_system()==0); }	

	GPUAPI void advance(){
		double H = min( max_timestep ,  _params.time_step );
		double pos = 0, vel = 0;

		if( is_in_body_component_grid() )
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

			p_i = p_im2 + 2.0 * h * v_im1;
			v_i = v_im2 + 2.0 * h * a_im1;
		}
		double a_i = calcForces.acc(ij,b,c,p_i,v_i);

		pos = 0.5 * ( p_i + p_im1 + h * v_i );
		vel = 0.5 * ( v_i + v_im1 + h * a_i );


		// Finalize the step
		if( is_in_body_component_grid() )
			sys[b][c].pos() = pos , sys[b][c].vel() = vel;
		if( is_first_thread_in_system() ) 
			sys.time() += H;
	}
};


} } } // End namespace bppt :: gpu :: swarm
