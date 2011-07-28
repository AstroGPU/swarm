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

#include "bppt.hpp"
#include "helpers.hpp"
#include "gravitation.hpp"


namespace swarm {

namespace gpu {
namespace bppt {

template< template<class T> class Propagator, class _Stopper >
class generic: public integrator {
	typedef integrator base;
	typedef  _Stopper stopper_t;
	typedef  typename Propagator<params_t<3> >::params prop_params_t;
	private:
	double _time_step;
	int _iteration_count;
	stopper_t _stopper;
	prop_params_t _prop_params;

	public:
	generic(const config& cfg): base(cfg),_time_step(0.001), _stopper(cfg),_prop_params(cfg) {
		if(!cfg.count("time step")) ERROR("Integrator gpu_generic requires a timestep ('time step' keyword in the config file).");
		_time_step = atof(cfg.at("time step").c_str());
	}

	virtual void launch_integrator() {
		_iteration_count = _destination_time / _time_step;
		launch_templatized_integrator(this);
	}


	template<class T>
	__device__ void kernel(T a){
		if(sysid()>=_dens.nsys()) return;

		// References to Ensemble and Shared Memory
		ensemble::SystemRef sys = _dens[sysid()];
		typedef typename Gravitation<T::n>::shared_data grav_t;
		Gravitation<T::n> calcForces(sys,*( (grav_t*) system_shared_data_pointer(a) ) );

		// Local variables
		const int nbod = T::n;
		// Body number
		int b = thread_body_idx(nbod);
		// Component number
		int c = thread_component_idx(nbod);
		int ij = thread_in_system();
		bool body_component_grid = (b < nbod) && (c < 3);
		bool first_thread_in_system = thread_in_system() == 0;


		// local variables
		typename stopper_t::tester stopper_tester = _stopper.get_tester(sys,*_log) ;
		Propagator<T> prop(_prop_params,sys,calcForces);
		prop.b = b;
		prop.c = c;
		prop.ij = ij;
		prop.body_component_grid = body_component_grid;
		prop.first_thread_in_system = first_thread_in_system;

		////////// INTEGRATION //////////////////////
		//
		prop.init();


		for(int iter = 0 ; (iter < _iteration_count) && sys.active() ; iter ++ ) {

			prop.advance();

			if( first_thread_in_system ) 
				sys.active() = ! stopper_tester() ;

			__syncthreads();
		}

	}


};


}
}
}
