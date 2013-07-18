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

/*! \file hermite_integrator.hpp
 *   \brief Defines and implements \ref swarm::gpu::bppt::hermite class - the 
 *          GPU implementation of PEC2 Hermite integrator.
 *
 */

#include "swarm/common.hpp"
#include "swarm/gpu/bppt.hpp"

namespace swarm { namespace gpu { namespace bppt {

/*! GPU implementation of PEC2 Hermite integrator
 * \ingroup integrators
 *
 */
template< class Monitor , template<class T> class Gravitation >
class hermite: public integrator {
	typedef integrator base;
	typedef Monitor monitor_t;
	typedef typename monitor_t::params mon_params_t;
	private:
	double _time_step;
	mon_params_t _mon_params;

public: //! Construct for class hermite integrator
	hermite(const config& cfg): base(cfg),_time_step(0.001), _mon_params(cfg) {
		_time_step =  cfg.require("time_step", 0.0);
	}

	virtual void launch_integrator() {
		launch_templatized_integrator(this);
	}
	
	
        //! Define the number of thread per system
	template<class T>
	static GENERIC int thread_per_system(T compile_time_param){
		const int grav = Gravitation<T>::thread_per_system();
		const int moni = Monitor::thread_per_system(compile_time_param);
		return max( grav, moni);
	}

        //! Define the amount of shared memory per system
	template<class T>
	static GENERIC int shmem_per_system(T compile_time_param){
		const int grav = Gravitation<T>::shmem_per_system();
		const int moni = Monitor::shmem_per_system(compile_time_param);
		return max( grav, moni);
	}	

        //! Convert internal coordinates to std coordinates
        GPUAPI void convert_internal_to_std_coord() {} 
        //! Convert std coordinates to internal coordinates
        GPUAPI void convert_std_to_internal_coord() {}  

	template<class T>
	__device__ void kernel(T compile_time_param){

		if(sysid()>=_dens.nsys()) return;
		
		// References to Ensemble and Shared Memory
		typedef Gravitation<T> Grav;
		ensemble::SystemRef sys = _dens[sysid()];
		
		typedef typename Monitor::shared_data data_t;
		monitor_t montest(_mon_params,sys,*_log, *((data_t*) system_shared_data_pointer(this,compile_time_param))) ;
				
		typedef typename Grav::shared_data grav_t;
		Grav calcForces(sys,*( (grav_t*) system_shared_data_pointer(this,compile_time_param) ) );
		
		

		// Local variables
		const int nbod = T::n;
		// Body number
		const int b = thread_body_idx(nbod);
		// Component number
		const int c = thread_component_idx(nbod);

		// local variables
		montest.init( thread_in_system() );
		__syncthreads();

		// local information per component per body
		double pos = 0.0, vel = 0.0 , acc0 = 0.0, jerk0 = 0.0;
		if( (b < nbod) && (c < 3) )
			{ pos = sys[b][c].pos(); vel = sys[b][c].vel(); }


		////////// INTEGRATION //////////////////////

		/// Calculate acceleration and jerk
		calcForces(thread_in_system(),b,c,pos,vel,acc0,jerk0);

		for(int iter = 0 ; (iter < _max_iterations) && sys.is_active() ; iter ++ ) 
		{
			double h = _time_step;

			if( sys.time() + h > _destination_time ) {
				h = _destination_time - sys.time();
			}

			
			/// Initial Evaluation
			/// Calculate the forces
			//calcForces(thread_in_system(),b,c,pos,vel,acc0,jerk0);

			// Predict 
			pos = pos +  h*(vel+(h*0.5)*(acc0+(h/3.0)*jerk0));
			vel = vel +  h*(acc0+(h*0.5)*jerk0);

			double pre_pos = pos, pre_vel = vel;

			double acc1,jerk1;
			{
				// Evaluation
				calcForces(thread_in_system(),b,c,pos,vel,acc1,jerk1);
				
				// Correct
#if 0 // OLD
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
#if 0 // OLD
				pos = pre_pos + (0.1-0.25) * (acc0 - acc1) * h * h - 1.0/60.0 * ( 7.0 * jerk0 + 2.0 * jerk1 ) * h * h * h;
				vel = pre_vel + ( -0.5 ) * (acc0 - acc1 ) * h -  1.0/12.0 * ( 5.0 * jerk0 + jerk1 ) * h * h;
#else
				pos = pre_pos + ((0.1-0.25) * (acc0 - acc1) - 1.0/60.0 * ( 7.0 * jerk0 + 2.0 * jerk1 ) * h )* h * h ;
				vel = pre_vel + (( -0.5 ) * (acc0 - acc1 ) -  1.0/12.0 * ( 5.0 * jerk0 + jerk1 ) * h ) * h ;
#endif
			}
			acc0 = acc1, jerk0 = jerk1;
			
			__syncthreads();
			montest.storeCurrentStat(thread_in_system());
			__syncthreads();
			
			/// Finalize the step
			if( (b < nbod) && (c < 3) )
				{ sys[b][c].pos() = pos; sys[b][c].vel() = vel; }
			if( thread_in_system()==0 ) 
				sys.time() += h;
			
			__syncthreads();
			montest( thread_in_system() );  
			__syncthreads();
			montest.init( thread_in_system() );
			__syncthreads();
			
			if( sys.is_active() && thread_in_system()==0 )  {
			    if( sys.time() >= _destination_time ) 
			    {	sys.set_inactive(); }
			}

			__syncthreads();


		}

	}


};

} } } // end namespace bppt :: integrators :: swarm
