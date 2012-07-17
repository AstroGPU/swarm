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
#include "../common.hpp"
#include "bppt.hpp"
#include "device_functions.h"

#define ASSUME_PROPAGATOR_USES_STD_COORDINATES 0

namespace swarm { namespace gpu { namespace bppt {

template<class T>
GENERIC const T& max3(const T& a, const T& b, const T& c){
	if( b > a )
		return c > b ? c : b;
	else
		return c > a ? c : a;
}

/**
 * \brief Generic integrator for rapid creation of new integrators.
 * \ingroup integrators
 *
 * The common functionality for a body-pair-per-thread GPU integrator is
 * put in this class. The core functionality of integration is provided by Propagator
 * template and the core functionality for logging and deactivation is provided by Monitor.
 * Inside the integration loop, Propagator and Monitor are used for the action work.
 *
 * Both Propagator and Monitor are supposed to have a struct "params". The params will be
 * initialized with the config. During integration, params is passed to the Propagator or
 * Monitor to set up initial values.
 *
 * Propagator is expected to expose three methods:
 *  - init: called before the integration loop (Usually empty)
 *  - advance: called inside the integration loop (Supposed to advance position, velocity and time)
 *  - shutdown: called after the integration loop (Usually empty)
 *
 * Monitor is expected to implement the "operator()" and act like a function object.
 * It is called inside the integration loop and is supposed to examine the state of the system.
 * Write notes to the log if necessary and deactivate the system if it does not need to integrated
 * anymore.
 *
 * For an example of how to use this Generic integrator refer to Euler integrator.
 *
 * \todo _time_step does not belog to here and it is not used here. The better place for 
 * _time_step is inside the Propagator
 *
 * An extention of this class proposed to take a Gravitation class as a template. It 
 * Is not implemented yet.
 *
 * one suggestion for the template line is:
 *  template< template<class T, class G> class Propagator, template<class L> class Monitor, class G >
 * and G is supposed to be the Gravitation class.
 */
template< template<class T,class G> class Propagator, class Monitor
	, template<class T> class Gravitation>
class generic: public integrator {
	typedef integrator base;

	//! Monitor is instantiated to write to GPU log.
	typedef Monitor monitor_t;

	//! Parameters of the monitor, should be initialized from config file
	typedef typename monitor_t::params mon_params_t;

	//! We don't really know number of bodies right now and it does not matter.
	//! parameters of the Propagator should be initialized with config file. But
	//! The only way to access Propagator class is to instantiate it with something.
	//
	typedef compile_time_params_t<3> defpar_t;
	typedef  typename Propagator< defpar_t, Gravitation<defpar_t> >::params prop_params_t;

	private:
	double _time_step;
	mon_params_t _mon_params;
	prop_params_t _prop_params;

	public:
	/**
	 * The integrator is initialized from the configuration
	 * The generic integrator does not require any configuration parameters.
	 * The configuration parameters are passed to Monitor::params and Propagator::params.
	 */
	generic(const config& cfg): base(cfg),_time_step(0.001), _mon_params(cfg),_prop_params(cfg) {
		_time_step = cfg.require("time_step", 0.0);
	}

	virtual void launch_integrator() {
		launch_templatized_integrator(this);
	}



	template<class T>
	static const int thread_per_system(T compile_time_param){
		const int grav = Gravitation<T>::thread_per_system();
		const int prop = Propagator<T,Gravitation<T> >::thread_per_system();
		const int moni = Monitor::thread_per_system(compile_time_param);
		return max3( grav, prop, moni);
	}

	template<class T>
	static GENERIC const int shmem_per_system(T compile_time_param){
		const int grav = Gravitation<T>::shmem_per_system();
		const int prop = Propagator<T,Gravitation<T> >::shmem_per_system();
		const int moni = Monitor::shmem_per_system(compile_time_param);
		return max3( grav, prop, moni);
	}

  //         __device__ void convert_internal_to_std_coord() {} ;
  //         __device__ void convert_std_to_internal_coord() {};

	/**
	 * \brief Integrator the system using the provided Propagator and Monitor.
	 *  This is the meat of things. We put all the common functionality in here
	 *  It basically initializes shared memory. Creates a gravitation object.
	 *  Finds the system to integrate. Runs the main integration loop.
	 *  Inside the loop, Propagator and Monitor are called alternatively
	 *
	 */
	template<class T>
	__device__ void kernel(T compile_time_param){
		if(sysid()>=_dens.nsys()) return;

		typedef Gravitation<T> GravitationInstance;

		// References to Ensemble and Shared Memory
		ensemble::SystemRef sys = _dens[sysid()];
		typedef typename GravitationInstance::shared_data grav_t;
		GravitationInstance calcForces(sys,*( (grav_t*) system_shared_data_pointer(this,compile_time_param) ) );

		/////////// Local variables /////////////
		const int nbod = T::n;               // Number of Bodies
		int b = thread_body_idx(nbod);       // Body id
		int c = thread_component_idx(nbod);  // Component id (x=0,y=1,z=2)
		int ij = thread_in_system();         // Pair id

		// Thread barrier predicates
		//		bool body_component_grid = (b < nbod) && (c < 3);         // Barrier to act on bodies and components
		//		bool first_thread_in_system = (thread_in_system() == 0);  // Barrier to select only the first thread


		// Setting up Monitor
		monitor_t montest(_mon_params,sys,*_log) ;

		// Setting up Propagator
		Propagator<T,GravitationInstance> prop(_prop_params,sys,calcForces);
		prop.b = b;
		prop.c = c;
		prop.ij = ij;
		//		prop.body_component_grid = body_component_grid;
		//		prop.first_thread_in_system = first_thread_in_system;

		////////// INTEGRATION //////////////////////
		
		prop.init();
		__syncthreads();


		for(int iter = 0 ; (iter < _max_iterations) && sys.is_active() ; iter ++ ) {

			prop.max_timestep = _destination_time - sys.time();
			prop.advance();
			__syncthreads();

			bool thread_needs_std_coord  = false;
			bool using_std_coord = false;


#if ASSUME_PROPAGATOR_USES_STD_COORDINATES
			montest( thread_in_system() );
#else
			thread_needs_std_coord = montest.pass_one( thread_in_system() );
#if (__CUDA_ARCH__ >= 200) 
			// requires arch=compute_sm_20
			bool block_needs_std_coord = syncthreads_or((int)(thread_needs_std_coord));
#else
#warning  Need to make this work for pre-Fermi GPUs.  For now just setting true!
			//			void *ptr_shared_void = calcForces.unused_shared_data_pointer(system_per_block_gpu());
			//			int *ptr_shared_int = static_cast<int *>(ptr_shared_void);
			//			//			int put_back = *ptr_shared_int;
			//						*ptr_shared_int = 0;
			//			//			atomicOr(ptr_shared_int,thread_needs_std_coord);
			//			  //			if(thread_needs_std_coord) (*ptr_shared_int)++;
			//			__syncthreads();
			//			bool block_needs_std_coord = static_cast<bool>(*ptr_shared_int);
			//			*ptr_shared_int = put_back;
			bool block_needs_std_coord = true;
#endif			
			if(block_needs_std_coord) 
			  { 
			    prop.convert_internal_to_std_coord(); 
			    using_std_coord = true; 
			  }

			__syncthreads();
			int new_state = montest.pass_two ( thread_in_system() );

			if( montest.need_to_log_system() && (thread_in_system()==0) )
			  { log::system(*_log, sys); }
			    
			__syncthreads();
			if(using_std_coord)
			  {
			    prop.convert_std_to_internal_coord();
			    using_std_coord = false;
			  }
#endif
			__syncthreads();			  

			if( sys.is_active() && prop.is_first_thread_in_system() )
			  {
			    if( sys.time() >= _destination_time ) 
			      { sys.set_inactive();     }
			  }
			__syncthreads();

		}

		prop.shutdown();

	}


};


} } } // End namespaces bppt, gpu, swarm
