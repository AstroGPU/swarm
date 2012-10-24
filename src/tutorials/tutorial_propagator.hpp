/*************************************************************************
 * Copyright (C) 2009-2010 by Eric Ford & the Swarm-NG Development Team  *
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

/** \file tutorial_propagator.hpp
 *  \brief Defines \ref TutorialPropagator class --
 *         implements a tutorial for making a propagator.
 *
 */

/*
 *  This is a simple tutorial used in doxygen pages
 *  should go through program2doxygen before it
 *  can be used by doxygen.
 *  
 *
 */
// \page TutorialPropagator Tutorial for Making a Propagator
//
// A propagator class implements a device function that advance 
// one system by one time step for one system (or at least attempts one 
// timestep). These can be readiliy combined with the generic_integrator 
// to quickly provide a new GPU-based integration algorithm.
//
// This is the header file where you define the new propagator
//
#include "swarm/swarmplugin.h"


struct TutorialPropagatorParams {
	double time_step;
	TutorialPropagatorParams(const config& cfg){
		time_step = cfg.require("time_step", 0.0);
	}
};

template<class T,class Gravitation>
struct TutorialPropagator {
	typedef TutorialPropagatorParams params;
	const static int nbod = T::n;

	params _params;


	// Runtime variables
	ensemble::SystemRef& sys;
	Gravitation& calcForces;
	int b;
	int c;
	int ij;
	bool body_component_grid;
	bool first_thread_in_system;
	double max_timestep;


	GPUAPI TutorialPropagator(const params& p,ensemble::SystemRef& s,
			Gravitation& calc)
		:_params(p),sys(s),calcForces(calc){}

	GPUAPI void init()  { }

	GPUAPI void shutdown() { }

        GPUAPI void convert_internal_to_std_coord() {} 
        GPUAPI void convert_std_to_internal_coord() {}

	static GENERIC int thread_per_system(){
		return nbod * 3;
	}

	static GENERIC int shmem_per_system() {
		 return 0;
	}

	GPUAPI bool is_in_body_component_grid()
        { return  ((b < nbod) && (c < 3)); }	

	GPUAPI bool is_in_body_component_grid_no_star()
        { return ( (b!=0) && (b < nbod) && (c < 3) ); }	

	GPUAPI bool is_first_thread_in_system()
        { return (thread_in_system()==0); }	

	GPUAPI void advance(){
		double pos = 0.0, vel = 0.0;
		double acc = 0.0, jerk = 0.0;
		
		if( is_in_body_component_grid() )
			pos = sys[b][c].pos() , vel = sys[b][c].vel();


		// First step of integration: calculate the accelartion 
		// (second derivative of position) and jerk (third derivative
		// of position). The gravitation
		// algorithm needs all the specified parameters. The acceleration and
		// jerk are returned in the last two variables (acc,jerk).
		calcForces(ij,b,c,pos,vel,acc,jerk);

		// For more precise integrations, we want to end exactly at
		// the destination_time, that means that h cannot be
		// greater than max_timestep that is allowed by the wrapper.
		double h = min(_params.time_step, max_timestep);

		// For this simple integrator, we use explicit Euler integration
		// equations. More complex equations can be used in practice.
		pos = pos +  h*(vel+(h*0.5)*(acc+(h/3.0)*jerk));
		vel = vel +  h*(acc+(h*0.5)*jerk);


		// Finalize the step
		if( is_in_body_component_grid() )
			sys[b][c].pos() = pos , sys[b][c].vel() = vel;
		if( is_first_thread_in_system() ) 
			sys.time() += h;
	}
};
