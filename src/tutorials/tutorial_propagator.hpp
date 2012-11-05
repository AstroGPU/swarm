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
// timestep). These can be readily combined with the generic_integrator 
// to quickly provide a new GPU-based integration algorithm.
//
// This is the header file where you define the new propagator
//
#include "swarm/swarmplugin.h"

// We have to create a separate data structure that holds the parameters
// for the propagator. This data structure is initialized with the
// configuration object when the propagator plugin is loaded. 
// 
struct TutorialPropagatorParams {
	double time_step;
	TutorialPropagatorParams(const config& cfg){
		time_step = cfg.require("time_step", 0.0);
	}
};

// This is the actual class that represents our propagator, we cannot
// initialize variables in the class because this class is instantiated
// at the place when it needs to be used. The propagator class is
// parametrized by number of bodies (class T contains it) and an
// implementation of Gravitational force calculation algorithm.
template<class T,class Gravitation>
class TutorialPropagator {
	public:
	
	// This will give the Generic integrator an idea about the struct
	// type we used for the parameters
	typedef TutorialPropagatorParams params;
	
	// We get the number of bodies at the compile-time. The propagator
	// should use this value for number of bodies instead of getting
	// it from the ensemble.
	const static int nbod = T::n;

	// We define variables that we use throughout the integration
	// _params contains all the parameters from the configuration. 
	// sys is the system that we have to operate on. and calcForces
	// is the implementation of force calculation algorithm
	// 
	// The constructure initializez these member variables from what
	// is provided.
	private:
	params _params;
	ensemble::SystemRef& sys;
	Gravitation& calcForces;
	
	GPUAPI TutorialPropagator(const params& p,ensemble::SystemRef& s,
			Gravitation& calc)
		:_params(p),sys(s),calcForces(calc){}

    // These are the variables that are generally used in the
    // integrators, we receive these variable in this way from
    // the integrator.
    // b, c and ij are define our work based on the thread id. b is
    // the number of body, c is the component number (0,1,2). and
    // ij is the pair number that is passed to calcForces
	public:
	int b;
	int c;
	int ij;
	
	// body_component_grid and first_thread_in_system are useful
	// predicates when we want to exclude some threads from updating
	// the data structures.
	bool body_component_grid;
	bool first_thread_in_system;
	
	// max_timestep is set by the generic_integrator, this is the biggest
	// time step that we are allowed to take. Usually it is only bound
	// by the destination_time and the default implementation uses
	// destination_time - sys.time().(but it can be used otherwise).
	double max_timestep;



	// init function is executed before entering the integration loop
	// a propagator can set-up data structure if needed.
	// shutdown is executed right after the integration loop. So it
	// do the clean-up.
	GPUAPI void init()  { }
	GPUAPI void shutdown() { }

	// These functions are only used if the propagator uses a coordinate
	// system other than the default.
	GPUAPI void convert_internal_to_std_coord() {} 
	GPUAPI void convert_std_to_internal_coord() {}

	// propagator can use arbitrary number of systems and may use
	// some shared memory, but it should be reported here so the 
	// launcher can initialize it.
	static GENERIC int thread_per_system(){ return nbod * 3; }
	static GENERIC int shmem_per_system() { return 0;        }



	// The advance function is called within the integration loop
	// The main purpose of the advance function is to integrate 
	// the system and advance the system in time. 
	//
	// The usual implemnation consist of sampling accleration (and jerk
	// if needed) at one or more points around the current system time
	// and extrapolate position and velocities.
	GPUAPI void advance(){
		// we define the local values just for more readable code.
		double pos = 0.0, vel = 0.0;
		double acc = 0.0, jerk = 0.0;
		
		// we have to use the predicate so we do not go out of bounds 
		// of the array.
		if( body_component_grid ) pos = sys[b][c].pos() , vel = sys[b][c].vel();


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
		pos = pos +  h * ( vel + (h*0.5) * (acc + (h/3.0)*jerk ) );
		vel = vel +  h * ( acc + (h*0.5) * jerk );


		// Finalize the step: save the position and velocities back to the ensemble
		// data structure and advance the system time.
		if(   body_component_grid  )  sys[b][c].pos() = pos , sys[b][c].vel() = vel;
		if( first_thread_in_system )  sys.time() += h;
	}
};

