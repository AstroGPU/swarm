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

/*! \file hermite_adap.hpp
 *   \brief Defines and implements \ref swarm::gpu::bppt::hermite_adap class - the  
 *          GPU implementation of PEC2 Hermite integrator with adaptive time step.
 *
 */

#include "swarm/common.hpp"
#include "swarm/gpu/bppt.hpp"

namespace swarm { namespace gpu { namespace bppt {

/*! GPU implementation of PEC2 Hermite integrator w/ adaptive time step
 * \ingroup integrators
 *
 */
template< class Monitor , template<class T> class Gravitation >
class hermite_adap: public integrator {
	typedef integrator base;
	typedef Monitor monitor_t;
	typedef typename monitor_t::params mon_params_t;
	private:
	double _time_step_factor, _min_time_step;
	mon_params_t _mon_params;
 
public:  //! constructor for hermite_adap
        hermite_adap(const config& cfg): base(cfg),_time_step_factor(0.001),_min_time_step(0.001), _mon_params(cfg) {         
		_time_step_factor =  cfg.require("time_step_factor", 0.0);
		_min_time_step =  cfg.require("min_time_step", 0.0);
	}

        template<class T> //! Shared memory size
	static GENERIC int shmem_per_system(T compile_time_param){
		return sizeof(SystemSharedData<T>)/SHMEM_CHUNK_SIZE;
	}

	virtual void launch_integrator() {
		launch_templatized_integrator(this);
	}

	template<class T> //! Date structure for system shared data
	struct SystemSharedData {   
  		typedef Gravitation<T> Grav;
		typename Grav::shared_data gravitation;
		DoubleCoalescedStruct<SHMEM_CHUNK_SIZE> time_step_factor[T::n];
	};

        GPUAPI void convert_internal_to_std_coord() {} 
        GPUAPI void convert_std_to_internal_coord() {}  

  template<class T>  //! calculate adaptive time steps
	__device__ double calc_adaptive_time_step(T compile_time_param, SystemSharedData<T>& shared, const double acc, const double jerk)
		{
		// Body number
		int b = thread_body_idx(T::n);
		// Component number
		int c = thread_component_idx(T::n);

		//! Put accelerations and jerks for each body and component into shared memory
		if( (b < T::n) && (c < 3) ) {
		    shared.gravitation[b][c].acc() = acc*acc;
		    shared.gravitation[b][c].jerk() = jerk*jerk;
	    }
		__syncthreads();
		//! calculate sum of squares of each component for each body
		//! store ratio in shared memory
		if( (b < T::n) && (c==0) ) {
		    double acc_mag_sq = shared.gravitation[b][0].acc()+shared.gravitation[b][1].acc()+shared.gravitation[b][2].acc();
		    double jerk_mag_sq = shared.gravitation[b][0].jerk()+shared.gravitation[b][1].jerk()+shared.gravitation[b][2].jerk();
		    shared.time_step_factor[b].value() = jerk_mag_sq/acc_mag_sq;
		}
		__syncthreads();
		if( thread_in_system() == 0 ) {
			double tf = shared.time_step_factor[0].value();
		    for(int bb=1;bb<T::n;++bb)
		      	tf += shared.time_step_factor[bb].value();
		    shared.time_step_factor[0].value() = rsqrt(tf)*_time_step_factor+_min_time_step;
		}
		__syncthreads();
		return shared.time_step_factor[0].value();
	}		



	template<class T>
	__device__ void kernel(T compile_time_param){

		if(sysid()>=_dens.nsys()) return;
		//! References to Ensemble and Shared Memory
		typedef Gravitation<T> Grav;
		ensemble::SystemRef sys = _dens[sysid()];
		SystemSharedData<T>& shared_data = *(SystemSharedData<T>*) system_shared_data_pointer(this, compile_time_param);
		//! Calculate the forces
		Grav calcForces(sys, shared_data.gravitation );

		// Local variables
		const int nbod = T::n;
		// Body number
		const int b = thread_body_idx(nbod);
		// Component number
		const int c = thread_component_idx(nbod);
//		const bool body_component_grid = (b < nbod) && (c < 3);


		// local variables
		monitor_t montest(_mon_params,sys,*_log) ;

		// local information per component per body
		double pos = 0.0, vel = 0.0 , acc0 = 0.0, jerk0 = 0.0;
		if( (b < T::n) && (c < 3) )
			pos = sys[b][c].pos() , vel = sys[b][c].vel();

		montest( thread_in_system() );

		////////// INTEGRATION //////////////////////

		//! Calculate acceleration and jerk
		calcForces(thread_in_system(),b,c,pos,vel,acc0,jerk0);

		for(int iter = 0 ; (iter < _max_iterations) && sys.is_active() ; iter ++ ) {
			// Since h is the time step that is used for the step it makes more sense to
			// to calculate time step and assign it to h
			double h = calc_adaptive_time_step(compile_time_param, shared_data ,acc0,jerk0);

			if( sys.time() + h > _destination_time ) {
				h = _destination_time - sys.time();
			}

			
			/// Initial Evaluation, it can be omitted for faster computation
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
#if 0 // OLD
				pos = pre_pos + (0.1-0.25) * (acc0 - acc1) * h * h - 1.0/60.0 * ( 7.0 * jerk0 + 2.0 * jerk1 ) * h * h * h;
				vel = pre_vel + ( -0.5 ) * (acc0 - acc1 ) * h -  1.0/12.0 * ( 5.0 * jerk0 + jerk1 ) * h * h;
#endif
				pos = pre_pos + ( (0.1-0.25) * (acc0 - acc1) - 1.0/60.0 * ( 7.0 * jerk0 + 2.0 * jerk1 ) * h) * h * h;
				vel = pre_vel + (( -0.5 ) * (acc0 - acc1 ) -  1.0/12.0 * ( 5.0 * jerk0 + jerk1 ) * h )* h ;
			}
			{
				// Evaluation
				calcForces(thread_in_system(),b,c,pos,vel,acc1,jerk1);
				
				// Correct
#if 0 // OLD
				pos = pre_pos + (0.1-0.25) * (acc0 - acc1) * h * h - 1.0/60.0 * ( 7.0 * jerk0 + 2.0 * jerk1 ) * h * h * h;
				vel = pre_vel + ( -0.5 ) * (acc0 - acc1 ) * h -  1.0/12.0 * ( 5.0 * jerk0 + jerk1 ) * h * h;
#endif
				pos = pre_pos + ((0.1-0.25) * (acc0 - acc1) - 1.0/60.0 * ( 7.0 * jerk0 + 2.0 * jerk1 ) * h )* h * h ;
				vel = pre_vel + (( -0.5 ) * (acc0 - acc1 ) -  1.0/12.0 * ( 5.0 * jerk0 + jerk1 ) * h ) * h ;
			}
			acc0 = acc1, jerk0 = jerk1;

			/// Finalize the step
			if( (b < T::n) && (c < 3) )
				sys[b][c].pos() = pos , sys[b][c].vel() = vel;
			if( thread_in_system()==0 ) 
				sys.time() += h;

			montest( thread_in_system() );

			if( sys.is_active() && thread_in_system()==0  )  {
			   if( sys.time() >= _destination_time ) 
			    {  	sys.set_inactive();    }

					}
			__syncthreads();


		}

	}


};


} } } // end namespace bppt :: integrators :: swarm
