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

/** \file tutorial_integrator.hpp
 *  \brief Defines \ref TutorialIntegrator class - 
 *         implements a tutorial for making a integrator.
 *
 */

// @page TutorialIntegrator Tutorial for implementing an Integrator
// If you are planning to implement a new ODE integration algorithm for
// use with Swarm and a propagator is not enough for your needs, you
// may consider writing an integrator. Most of the code consists of 
// structures that other components of Swarm expect from an integrator.
//
// Actual mathematical equations embedded in the kernel function.
//
// First you need to include the headers for the parent classes.
#include "swarm/common.hpp"
#include "swarm/gpu/bppt.hpp"

// We need to use classes in Swarm namespace as well as classes in
// Swarm body-pair-per-thread namespace.
using swarm::config;
using swarm::ensemble;
using namespace swarm::gpu::bppt;

// The integrator does not have to be a template, but defining it as
// a template makes it easier to use different Monitors and different
// Gravitational force calculation algorithms
template< class Monitor , template<class T> class Gravitation >
// Name the integrator, and implement swarm::gpu::bppt::integrator
class TutorialIntegrator: public integrator {
	
	// Some convenience aliases, just to follow the conventions
	typedef integrator base;
	typedef Monitor monitor_t;
	typedef typename monitor_t::params mon_params_t;
	
	// Configuration paramaters that are loaded during when the integrator
	// is created. Values for these parameters are usually specified
	// in a config file.
	private:
	double _time_step;
	mon_params_t _mon_params;

	// The configurations are loaded when the integrator is first created.
	// The integrator should also configure the monitor by initializing the
	// _mon_params and passing the cfg parameter to it.
	public:
	TutorialIntegrator(const config& cfg): base(cfg),_time_step(0.001), _mon_params(cfg) {
		_time_step =  cfg.require("time_step", 0.0);
	}

	// Every integrator should implement the launch_integrator method. 
	// It is an abstract method, C++ will require you to implement it.
	// for most integrators that are templatized by the number of bodies.
	// it is easire to call launch_templatized_integrator on self. 
	// launch_templatized integrator will call a version of the
	// kernel function (see below) optimized for the specific number of
	// bodies.
	virtual void launch_integrator() {
		launch_templatized_integrator(this);
	}

		// These two function are used by some monitors. If you use
		// a coordinate system other than Cartesian, you may need
		// to convert it to Cartesian and back when these functions are called.
        GPUAPI void convert_internal_to_std_coord() {} 
        GPUAPI void convert_std_to_internal_coord() {}

	// The CUDA kernel that contains our algorithm is defined in this
	// `kernel` function. The template contains the compile-time parameters
	// which this kernel is being optimized for. Usaually it only contains
	// the number of bodies in a system.
	//
	// There are no dynamic parameters, since everything else (ensemble
	// and configuration parameters) is already
	// specified in the member variables of current class.
	template<class T>
	__device__ void kernel(T compile_time_param){

		// Usally there are more threads than the number of systems, so
		// we have to guard.
		if(sysid()>=_dens.nsys()) return;
		
		// We define an auxiliary variable for our system since our thread
		// only handles one system from the whole ensemble.
		ensemble::SystemRef sys = _dens[sysid()];

		// The Gravitation class needs to be optimized with the compile
		// time parameter. It also needs to use shared memory. So we 
		// use the helper function system_shared_data_pointer to retrieve
		// the shared memory area dedicated for our system.
		typedef Gravitation<T> Grav;
		typedef typename Grav::shared_data grav_t;
		Grav calcForces(sys,*( (grav_t*) system_shared_data_pointer(this,compile_time_param) ) );

		// We initialize the monitor and pass it the configuration
		// parameters and pointers to the current system and logging 
		// object.
		monitor_t montest(_mon_params,sys,*_log) ;

		// Here we use some variables to make our code look clearer.
		//
		// Number of bodies come from the compile time parameter since
		// our code is being optimized for that number of bodies.
		const int nbod = T::n;
		// This thread is assigned one component (x, y or z) of exactly
		// one body to work with, we use
		// the utility functions to find which body and which component
		// is assigned to current thread.
		const int b = thread_body_idx(nbod);
		const int c = thread_component_idx(nbod);


		// As a part of preparing to integrate, we load the values
		// to local variables. Local variables are usually represent
		// registers in native code. We always need to use 
		// the guards (b < nbod)&&(c < 3) because there are some threads
		// for which b and c are not valid. 
		double pos = 0.0, vel = 0.0 , acc = 0.0, jerk = 0.0;
		if( (b < nbod) && (c < 3) )
			{ pos = sys[b][c].pos(); vel = sys[b][c].vel(); }

		// We need to test with the monitor before we proceed since
		// the system may not get stopped before getting integrated
	    montest( thread_in_system() );


		// We use _max_iteration to avoid falling into infinite loop
		// (or very long CUDA calls). Otherwise, the loop only ends
		// when the system becomes inactive. (due to reaching the destination_time
		// or triggerring a stopping condition).
		for(int iter = 0 ; (iter < _max_iterations) && sys.is_active() ; iter ++ ) 
		{
			

			// First step of integration: calculate the accelartion 
			// (second derivative of position) and jerk (third derivative
			// of position). The gravitation
			// algorithm needs all the specified parameters. The acceleration and
			// jerk are returned in the last two variables (acc,jerk).
			calcForces(thread_in_system(),b,c,pos,vel,acc,jerk);

			// For more precise integrations, we want to end exactly at
			// the destination_time, that means that h cannot be
			// greater than _desitantion_time - current_system_time
			double h = min(_destination_time - sys.time(), _time_step);
				 
			// For this simple integrator, we use explicit Euler integration
			// equations. More complex equations can be used in practice.
			pos = pos +  h*(vel+(h*0.5)*(acc+(h/3.0)*jerk));
			vel = vel +  h*(acc+(h*0.5)*jerk);

			// Write the positions and velocities back to memory
			// for the monitor testing to use it
			if( (b < nbod) && (c < 3) )
				{ sys[b][c].pos() = pos; sys[b][c].vel() = vel; }
			
			// Only the first thread should advance the current 
			// system time
			if( thread_in_system()==0 ) 
				sys.time() += h;
				
			// We need to make sure that all the memory writes has
			// been done before we can proceed.
			__syncthreads();
			
			montest( thread_in_system() );  
			__syncthreads();
			
			// We also have to check that if we reached the destination_time
			// if that is the case, we deactivate system so we can 
			// get out of the loop.
			if( sys.is_active() && thread_in_system()==0 )  {
			    if( sys.time() >= _destination_time ) 
			    {	sys.set_inactive(); }
			}

			// We need to do resynchronize so all threads can see 
			// the changes in time and state of the system.
			__syncthreads();

		} // end of for loop
	} // end of kernel function
}; // end of class TutorialIntegrator

