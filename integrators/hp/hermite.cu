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

#include "hp.hpp"
#include "swarmlog.h"


namespace swarm {
namespace hp {

class hermite: public integrator {
	typedef integrator base;
	private:
	double _time_step;

	public:
	hermite(const config& cfg): base(cfg),_time_step(0.001) {
		if(!cfg.count("time step")) ERROR("Integrator gpu_hermite requires a timestep ('time step' keyword in the config file).");
		_time_step = atof(cfg.at("time step").c_str());
	}

	virtual void launch_integrator() {
		launch_templatized_integrator(this);
	}

	/*! Integrator Kernel to be run on GPU
	 *  
	 *
	 */
	 template<class T >
	__device__ void kernel(T a)  {

		const int nbod = T::n;

		/////////////// FETCH LOCAL VARIABLES ///////////////////

		int thr = thread_in_system();

		if(sysid() >= _gpu_ens->nsys()) { return; }
		ensemble::systemref sys ( (*_gpu_ens)[sysid()] );

		// Body/Component Grid
		// Body number
		int b = thr / 3 ;
		// Component number
		int c = thr % 3 ;
		bool body_component_grid = b < nbod;

		// i,j pairs Grid
		int ij = thr;

		// shared memory allocation
		extern __shared__ char shared_mem[];
		char*  system_shmem =( shared_mem + sysid_in_block() * integrator::shmem_per_system(nbod) );

		double t_start = sys.time(), t = t_start;
		// TODO: Change the way stopping is done
		double t_end = min(_destination_time,sys.time_end());
//		double t_end = min(t_start + _destination_time,sys.time_end());



		// local information per component per body
		double pos = 0, vel = 0 , acc0 = 0, jerk0 = 0;
		if( body_component_grid )
			pos = sys[b].p(c), vel = sys[b].v(c);


		////////// INTEGRATION //////////////////////

		// Calculate acceleration and jerk
		Gravitation<nbod> calcForces(sys,system_shmem);
		calcForces(ij,b,c,pos,vel,acc0,jerk0);

		unsigned int iter=0;     // Make sure don't get stuck in infinite loop
		bool skip_loop = ( t + _time_step > t_end );
		while(!skip_loop)
		   {
//			double h = min(_time_step, t_end - t);
			double h = _time_step;
			
			// Initial Evaluation
			// WARNING:  Reusing results from call at end of loop
			///calcForces(ij,b,c,pos,vel,acc0,jerk0);

			// Predict 
			pos = pos +  h*(vel+(h*0.5)*(acc0+(h/3.)*jerk0));
			vel = vel +  h*(acc0+(h*0.5)*jerk0);

			double pre_pos = pos, pre_vel = vel;

			double acc1,jerk1;
			{
				// Evaluation
				calcForces(ij,b,c,pos,vel,acc1,jerk1);

				// Correct
				pos = pre_pos + (.1-.25) * (acc0 - acc1) * h * h - 1/60.0 * ( 7 * jerk0 + 2 * jerk1 ) * h * h * h;
				vel = pre_vel + ( -.5 ) * (acc0 - acc1 ) * h -  1/12.0 * ( 5 * jerk0 + jerk1 ) * h * h;
			}
			{
				// Evaluation
				calcForces(ij,b,c,pos,vel,acc1,jerk1);

				// Correct
				pos = pre_pos + (.1-.25) * (acc0 - acc1) * h * h - 1/60.0 * ( 7 * jerk0 + 2 * jerk1 ) * h * h * h;
				vel = pre_vel + ( -.5 ) * (acc0 - acc1 ) * h -  1/12.0 * ( 5 * jerk0 + jerk1 ) * h * h;
			}
			acc0 = acc1, jerk0 = jerk1;

			// Finalize the step
			t += h;
                        if( (t+h>t_end) || (iter>=_max_itterations_per_kernel_call) ) 
			   skip_loop = true;

			// Need to write these, so visible to calcForces next itteration (and upon exit)
			if( body_component_grid )
				sys[b].p(c) = pos, sys[b].v(c) = vel;

                        // WARNING: If needs_output were to read pos and vel from _gpu_ens, then it would need to call syncthreads itself
                        const bool sys_needs_output = (thr==0) ? log::needs_output(*_gpu_ens, t, sysid()) : 0;

                        if(sys_needs_output) // implicit if (thr == 0) 
                           {
		  	   sys.set_time(t);  // Only write to time if going to output
                           }
                        __syncthreads();   // Since need other threads within block to have completed writing to sys[b].p/v(c) before calling output_system 
                        if(sys_needs_output) // implicity if (thr == 0) 
		 	    {
			    log::output_system(*_gpu_log, *_gpu_ens, t, sysid());
			    }
		}

		if(thr == 0) 
			sys.set_time(t);
	     	++iter;	
	}

};

/*!
 * \brief Factory to create double/single/mixed hermite gpu integrator based on precision
 *
 * @param[in] cfg configuration class
 *
 * @return        pointer to integrator cast to integrator*
 */
extern "C" integrator *create_hp_hermite(const config &cfg)
{
	return new hermite(cfg);
}

}
}

