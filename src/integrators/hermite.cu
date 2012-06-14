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

#include "swarm/common.hpp"
#include "swarm/gpu/bppt.hpp"
#include "monitors/composites.hpp"
#include "monitors/stop_on_ejection.hpp"
#include "monitors/log_time_interval.hpp"


namespace swarm { namespace gpu { namespace bppt {

/*! GPU implementation of PEC2 Hermite integrator
 * \ingroup integrators
 *
 */
template< class Monitor >
class hermite: public integrator {
	typedef integrator base;
	typedef Monitor monitor_t;
	typedef typename monitor_t::params mon_params_t;
	private:
	double _time_step;
	mon_params_t _mon_params;

	public:
	hermite(const config& cfg): base(cfg),_time_step(0.001), _mon_params(cfg) {
		_time_step =  cfg.require("time_step", 0.0);
	}

	virtual void launch_integrator() {
		launch_templatized_integrator(this);
	}


        GPUAPI void convert_internal_to_std_coord() {} 
        GPUAPI void convert_std_to_internal_coord() {}

	template<class T>
	__device__ void kernel(T compile_time_param){

		if(sysid()>=_dens.nsys()) return;
		// References to Ensemble and Shared Memory
		ensemble::SystemRef sys = _dens[sysid()];
		typedef typename Gravitation<T::n>::shared_data grav_t;
		Gravitation<T::n> calcForces(sys,*( (grav_t*) system_shared_data_pointer(this,compile_time_param) ) );

		// Local variables
		const int nbod = T::n;
		// Body number
		const int b = thread_body_idx(nbod);
		// Component number
		const int c = thread_component_idx(nbod);

		// local variables
		monitor_t montest(_mon_params,sys,*_log) ;


		// local information per component per body
		double pos = 0.0, vel = 0.0 , acc0 = 0.0, jerk0 = 0.0;
		if( (b < nbod) && (c < 3) )
			{ pos = sys[b][c].pos(); vel = sys[b][c].vel(); }


//		if( thread_in_system()==0  )  {
		    montest( thread_in_system() );
//		    }

		////////// INTEGRATION //////////////////////

		// Calculate acceleration and jerk
		calcForces(thread_in_system(),b,c,pos,vel,acc0,jerk0);

		for(int iter = 0 ; (iter < _max_iterations) && sys.is_active() ; iter ++ ) 
		{
			double h = _time_step;

			if( sys.time() + h > _destination_time ) {
				h = _destination_time - sys.time();
			}

			
			// Initial Evaluation
			///calcForces(thread_in_system(),b,c,pos,vel,acc0,jerk0);

			// Predict 
			pos = pos +  h*(vel+(h*0.5)*(acc0+(h/3.0)*jerk0));
			vel = vel +  h*(acc0+(h*0.5)*jerk0);

			double pre_pos = pos, pre_vel = vel;

			double acc1,jerk1;
			{
				// Evaluation
				calcForces(thread_in_system(),b,c,pos,vel,acc1,jerk1);
				
				// Correct
#if 1 // OLD
				pos = pre_pos + (0.1-0.25) * (acc0 - acc1) * h * h - 1.0/60.0 * ( 7.0 * jerk0 + 2.0 * jerk1 ) * h * h * h;
				vel = pre_vel + ( -0.5 ) * (acc0 - acc1 ) * h -  1.0/12.0 * ( 5.0 * jerk0 + jerk1 ) * h * h;
#else
				pos = pre_pos + ( (0.1-0.25) * (acc0 - acc1) - 1.0/60.0 * ( 7.0 * jerk0 + 2.0 * jerk1 ) * h) * h * h;
				vel = pre_vel + (( -0.5 ) * (acc0 - acc1 ) -  1.0/12.0 * ( 5.0 * jerk0 + jerk1 ) * h )* h ;
#endif
			}
			{
				// Evaluation
				calcForces(thread_in_system(),b,c,pos,vel,acc1,jerk1);
				
				// Correct
#if 1 // OLD
				pos = pre_pos + (0.1-0.25) * (acc0 - acc1) * h * h - 1.0/60.0 * ( 7.0 * jerk0 + 2.0 * jerk1 ) * h * h * h;
				vel = pre_vel + ( -0.5 ) * (acc0 - acc1 ) * h -  1.0/12.0 * ( 5.0 * jerk0 + jerk1 ) * h * h;
#else
				pos = pre_pos + ((0.1-0.25) * (acc0 - acc1) - 1.0/60.0 * ( 7.0 * jerk0 + 2.0 * jerk1 ) * h )* h * h ;
				vel = pre_vel + (( -0.5 ) * (acc0 - acc1 ) -  1.0/12.0 * ( 5.0 * jerk0 + jerk1 ) * h ) * h ;
#endif
			}
			acc0 = acc1, jerk0 = jerk1;

			// Finalize the step
			if( (b < nbod) && (c < 3) )
				{ sys[b][c].pos() = pos; sys[b][c].vel() = vel; }
			if( thread_in_system()==0 ) 
				sys.time() += h;
			__syncthreads();
			montest( thread_in_system() );  
			__syncthreads();
			if( sys.is_active() && thread_in_system()==0 )  {
			    if( sys.time() >= _destination_time ) 
			    {	sys.set_inactive(); }
			}

			__syncthreads();


		}

	}


};


typedef gpulog::device_log L;
using namespace monitors;

integrator_plugin_initializer<hermite< stop_on_ejection<L> > >
	hermite_plugin("hermite");

integrator_plugin_initializer<hermite< stop_on_ejection_or_close_encounter<L> > >
	hermite_close_encounter_plugin("hermite_close_encounter");

integrator_plugin_initializer<hermite< log_time_interval<L> > >
	hermite_log_plugin("hermite_log");


} } } // end namespace bppt :: integrators :: swarm
