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

#include "hp/bppt.hpp"
#include "hp/helpers.hpp"
#include "hp/gravitation.hpp"


namespace swarm {
namespace hp {

namespace gpu {
namespace bppt {

class hermite: public integrator {
	typedef integrator base;
	private:
	double _time_step;
	int _iteration_count;

	public:
	hermite(const config& cfg): base(cfg),_time_step(0.001) {
		if(!cfg.count("time step")) ERROR("Integrator gpu_hermite requires a timestep ('time step' keyword in the config file).");
		_time_step = atof(cfg.at("time step").c_str());
	}

	virtual void launch_integrator() {
		_iteration_count = _destination_time / _time_step;
		launch_templatized_integrator(this);
	}


	template<class T>
	__device__ void kernel(T a){
		// References to Ensemble and Shared Memory
		ensemble::SystemRef sys = _dens[sysid()];
		extern __shared__ char shared_mem[];
		CoalescedStructArray< typename Gravitation<T::n>::shared_data
			, double, SHMEM_WARPSIZE > 
			grav_shared((typename Gravitation<T::n>::shared_data*) shared_mem,0);
		Gravitation<T::n> calcForces(sys,grav_shared[sysid_in_block()]);

		// Local variables
		const int nbod = T::n;
		// Body number
		int b = thread_body_idx(nbod);
		// Component number
		int c = thread_component_idx(nbod);
		int ij = thread_in_system();
		bool body_component_grid = (b < nbod) && (c < 3);

		double t = sys.time(); 
		double t_end = sys.time() + _destination_time;

		// local information per component per body
		double pos = 0, vel = 0 , acc0 = 0, jerk0 = 0;
		if( body_component_grid )
			pos = sys[b][c].pos() , vel = sys[b][c].vel();


		////////// INTEGRATION //////////////////////

		// Calculate acceleration and jerk
		calcForces(ij,b,c,pos,vel,acc0,jerk0);

		for(int iter = 0 ; iter < _iteration_count; iter ++ ) {
			double h = min(_time_step, t_end - t);

			
			// Initial Evaluation
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


		}
		sys.time() = t;
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
}
}
