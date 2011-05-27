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
#include "static_accjerk.hpp"
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
		double t_end = min(t_start + _destination_time,sys.time_end());

		// local information per component per body
		double pos = 0, vel = 0 , acc = 0, jerk = 0;
		if( body_component_grid )
			pos = sys[b].p(c), vel = sys[b].v(c);


		////////// INTEGRATION //////////////////////

		// Calculate acceleration and jerk
		Gravitation<nbod> calcForces(sys,system_shmem);
		calcForces(ij,b,c,pos,vel,acc,jerk);

		while(t < t_end){
			for(int k = 0; k < 2; k++)
			{
				double h = min(_time_step, t_end - t);
				double pos_old = pos, vel_old = vel, acc_old = acc,jerk_old = jerk;

				// Predict 
				pos = pos_old +  h*(vel_old+(h*0.5)*(acc+(h/3.)*jerk));
				vel = vel_old +  h*(acc+(h*0.5)*jerk);

				// Do evaluation and correction two times (PEC2)
				for(int l = 0; l < 2; l++)
				{

					// Calculate acceleration and jerk using shared memory
					calcForces(ij,b,c,pos,vel,acc,jerk);

					// Correct
					pos = pos_old + (h*0.5) * ( (vel_old + vel) 
							+ (h*7.0/30.)*( (acc_old-acc) + (h/7.) * (jerk_old+jerk)));
					vel = vel_old + (h*0.5) * ( (acc_old+acc) + (h/6.) * (jerk_old-jerk));

				}
				t += h;
			}

			if( body_component_grid )
				sys[b].p(c) = pos, sys[b].v(c) = vel;

			if(thr == 0) 
				if(log::needs_output(*_gpu_ens, t, sysid()))
				{
					sys.set_time(t);
					log::output_system(*_gpu_log, *_gpu_ens, t, sysid());
				}

		}

		if(thr == 0) 
			sys.set_time(t);

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

