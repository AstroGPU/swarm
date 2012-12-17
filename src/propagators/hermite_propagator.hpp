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

/*! \file hermite_propagator.hpp
 *   \brief Defines swarm::gpu::bppt::HermitePropagator - implements the GPU version 
 *          of hermite propagator. 
 *
 */


#include "swarm/swarmplugin.h"

namespace swarm {

namespace gpu {
namespace bppt {

/*! Paramaters for HermitePropagator
 * \ingroup propagator_parameters
 *
 */
struct HermitePropagatorParams {
	double time_step;
        //! Constructor
	HermitePropagatorParams(const config& cfg){
		time_step = cfg.require("time_step", 0.0);
	}
};

/*! GPU implementation of hermite propagator
 * It is of no practical use since @ref hermite integrator implements
 * the same functionaliy faster. It is only given for performance comparison
 * between equivalent propagators and integrators.
 * 
 * \ingroup propagators
 *
 */
template<class T,class Gravitation>
struct HermitePropagator {
	typedef HermitePropagatorParams params;
	const static int nbod = T::n;

	params _params;

	//! Runtime variables
	ensemble::SystemRef& sys;
	Gravitation& calcForces;
	int b;
	int c;
	int ij;
	double acc0, jerk0;
	bool body_component_grid;
	bool first_thread_in_system;
	double max_timestep;

        //! Constructor for HermitePropagator
	GPUAPI HermitePropagator(const params& p,ensemble::SystemRef& s,
			Gravitation& calc)
		:_params(p),sys(s),calcForces(calc){}

        //! Initialize the system
	GPUAPI void init()  {
		double pos,vel;
		if( is_in_body_component_grid() )
			pos = sys[b][c].pos() , vel = sys[b][c].vel();
		calcForces(ij,b,c,pos,vel,acc0,jerk0);
	}

	GPUAPI void shutdown() { }

        //! Conversation between the coordinate systems. 
        GPUAPI void convert_internal_to_std_coord() {} 
        GPUAPI void convert_std_to_internal_coord() {}

	static GENERIC int thread_per_system(){
		return nbod * 3;
	}

	static GENERIC int shmem_per_system() {
		 return 0;
	}

	__device__ bool is_in_body_component_grid()
        { return  ((b < nbod) && (c < 3)); }	

	__device__ bool is_in_body_component_grid_no_star()
        { return ( (b!=0) && (b < nbod) && (c < 3) ); }	

	__device__ bool is_first_thread_in_system()
        { return (thread_in_system()==0); }	

        //! Advance the time steps
	GPUAPI void advance(){
		double h = min(_params.time_step, max_timestep);
		double pos = 0.0, vel = 0.0;
		double acc = 0.0, jerk = 0.0;
		const double third = 1.0/3.0;
		
		if( is_in_body_component_grid() )
			pos = sys[b][c].pos() , vel = sys[b][c].vel();


		//! Predict 
		pos = pos +  h*(vel+(h*0.5)*(acc0+(h/3.0)*jerk0));
		vel = vel +  h*(acc0+(h*0.5)*jerk0);

		double pre_pos = pos, pre_vel = vel;

		double acc1,jerk1;
#pragma unroll
		for(int i = 0; i < 2 ; i++){
			//! Evaluation
			calcForces(ij,b,c,pos,vel,acc1,jerk1);
			
			//! Correct
			pos = pre_pos + ( (0.1-0.25) * (acc0 - acc1) - 1.0/60.0 * ( 7.0 * jerk0 + 2.0 * jerk1 ) * h) * h * h;
			vel = pre_vel + (( -0.5 ) * (acc0 - acc1 ) -  1.0/12.0 * ( 5.0 * jerk0 + jerk1 ) * h )* h ;
		}
		acc0 = acc1, jerk0 = jerk1;


		//! Finalize the step
		if( is_in_body_component_grid() )
			sys[b][c].pos() = pos , sys[b][c].vel() = vel;
		if( is_first_thread_in_system() ) 
			sys.time() += h;
	}
};

}
}
}

