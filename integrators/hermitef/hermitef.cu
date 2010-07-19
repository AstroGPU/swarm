/*************************************************************************
 * Copyright (C) 2010 by Swarm-NG Development Team                       *
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

/*! \file hermitef.cu
 *  \brief declares prop_hermitef for use with gpu_generic_integrator
*/

#include "swarm.h"
#include "hermitef.h"
#include "swarmlog.h"
#include "meta.hpp"


/// namespace for Swarm-NG library
namespace swarm {

/*!  
 *  \brief propagator class for Hermite integrator on GPU: Advance the system by one time step.
 *
 *  CPU state and interface. Will be instantiated on construction of gpu_generic_integrator object. 
 *  Keep any data that need to reside on the CPU here.
 */
struct prop_hermitef
{
	/*! 
	 * \brief GPU state and interface (per-grid). Will be passed as an argument to integration kernel. 
         *
	 * Any per-block read-only variables should be members of this structure.
         */
	struct gpu_t
	{
		//! per-block variables, time step
		double h;
	
		/*!
		 * \brief  GPU per-thread state and interface. Will be instantiated in the integration kernel. 
                 *
		 * Any per-thread variables should be members of this structure.
		 */
		struct thread_state_t
		{
			thread_state_t(const gpu_t &H, ensemble &ens, const int sys, double T, double Tend)
			{ }
		};



		template<int nbod>
		__device__ static void corrector(int i,const int& c,ensemble::systemref sysref
				,double (&pos)[3][nbod],double (&vel)[3][nbod],double (&acc)[3][nbod],double (&jerk)[3][nbod]
				,double (&acc_old)[3][nbod],double (&jerk_old)[3][nbod]
				,const double& h){
			pos[c][i] = sysref[i].p(c) + (h/2) * ( (sysref[i].v(c)+vel[c][i]) 
					+ (h*7.0/30)*( (acc_old[c][i]-acc[c][i]) + (h/7) * (jerk_old[c][i]+jerk[c][i])));
			vel[c][i] = sysref[i].v(c) + (h/2) * ( (acc_old[c][i]+acc[c][i]) + (h/6) * (jerk_old[c][i]-jerk[c][i]));
		}


		template<int nbod>
		__device__ static void predictor(int i,int c
				,double (&pos)[3][nbod],double (&vel)[3][nbod],double (&acc)[3][nbod],double (&jerk)[3][nbod]
				,double h){
			pos[c][i] +=  h*(vel[c][i]+(h/2)*(acc[c][i]+(h/3)*jerk[c][i]));
			vel[c][i] +=  h*(acc[c][i]+(h/2)*jerk[c][i]);
		}

		/**!
		 * helper function for accjerk_updater that operates on each component
		 * it gets scalar part of acceleration as input and calculates one component of
		 * acceleration and jerk at a time
         *
		 */
		template<int nbod>
		__device__ static void accjerk_updater_component(int i,int c
				,double dx[3],double dv[3],double scalar,double rv
				,double (&acc)[3][nbod],double (&jerk)[3][nbod]){
			acc[c][i] += dx[c]* scalar;
			jerk[c][i] += (dv[c] - dx[c] * rv) * scalar;

		}

		inline __device__ static double inner_product(const double a[3],const double b[3]){
			return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
		}

		inline __device__ static double sqr(const double &a){
			return a*a;
		}

		template<int nbod>
		__device__ static void accjerk_updater(int ij,ensemble::systemref sys
				,double (&pos)[3][nbod],double (&vel)[3][nbod],double (&acc)[3][nbod],double (&jerk)[3][nbod]){
			const int i = ij/nbod, j = ij%nbod;
			if(i < j){

				double dx[3] =  { pos[0][j]-pos[0][i],pos[1][j]-pos[1][i],pos[2][j]-pos[2][i]};
				double dv[3] =  { vel[0][j]-vel[0][i],vel[1][j]-vel[1][i],vel[2][j]-vel[2][i]};

				// computing scalar part of the acceleration
				double r2 =  dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2] ;
				double rv =  inner_product(dx,dv) * 3 / r2;
				double rinv = rsqrt(r2)  / r2;

				// vectorized part
				const double scalar_i = +rinv*sys[j].mass();
				accjerk_updater_component<nbod>(i,0,dx,dv,scalar_i,rv,acc,jerk);
				accjerk_updater_component<nbod>(i,1,dx,dv,scalar_i,rv,acc,jerk);
				accjerk_updater_component<nbod>(i,2,dx,dv,scalar_i,rv,acc,jerk);

				const double scalar_j = -rinv*sys[i].mass();
				accjerk_updater_component<nbod>(j,0,dx,dv,scalar_j,rv,acc,jerk);
				accjerk_updater_component<nbod>(j,1,dx,dv,scalar_j,rv,acc,jerk);
				accjerk_updater_component<nbod>(j,2,dx,dv,scalar_j,rv,acc,jerk);
			}
		}


		template<int nbod>
	    __device__ static void advance_internal(ensemble::systemref& sysref, double h){
			// Local storage that should be in registers
			double pos[3][nbod], vel[3][nbod], acc[3][nbod], jerk[3][nbod], acc_old[3][nbod], jerk_old[3][nbod];
			using namespace boost;

			if(sysref.nbod() == nbod) {

				// load position and velocity from memory
				for(int c = 0; c < 3 ; c++)
					#pragma unroll
					for(int i = 0; i < nbod;i++)
						pos[c][i] = sysref[i].p(c), vel[c][i] = sysref[i].v(c);
					

				for(int c = 0; c < 3 ; c++)
					#pragma unroll
					for(int i = 0; i < nbod;i++)
						acc_old[c][i] = 0, jerk_old[c][i] = 0;
				// Update Acc Jerk before predicting
				//Unroller<0,nbod*nbod>::step(bind<void>(&accjerk_updater<nbod>,boost::arg<1>(),sysref,pos,vel,acc_old,jerk_old));
				#pragma unroll
				for(int ij = nbod*nbod-1; ij >=0 ;ij--)
					accjerk_updater<nbod>(ij,sysref,pos,vel,acc_old,jerk_old);
				

				// Predict
				//Unroller<0,nbod>::step(bind<void>(predictor<nbod>,_1,0,pos,vel,acc,jerk,h));
				for(int c = 0; c < 3 ; c++)
					#pragma unroll
					for(int i = 0; i < nbod;i++)
						predictor<nbod>(i,c,pos,vel,acc_old,jerk_old,h);


				for(int c = 0; c < 3 ; c++)
					#pragma unroll
					for(int i = 0; i < nbod;i++)
						acc[c][i] = 0, jerk[c][i] = 0;
				// Update Acc Jerk
				//Unroller<0,nbod*nbod>::step(bind<void>(accjerk_updater<nbod>,_1,sysref,pos,vel,acc,jerk));
				#pragma unroll
				for(int ij = nbod*nbod-1; ij >=0 ;ij--)
					accjerk_updater<nbod>(ij,sysref,pos,vel,acc,jerk);

				// Correct
				//Unroller<0,nbod>::step(bind<void>(corrector<nbod>,_1,0,sysref,pos,vel,acc,jerk,acc_old,jerk_old,h));
				for(int c = 0; c < 3 ; c++)
					#pragma unroll
					for(int i = 0; i < nbod;i++)
						corrector<nbod>(i,c,sysref,pos,vel,acc,jerk,acc_old,jerk_old,h);

				for(int c = 0; c < 3 ; c++)
					#pragma unroll
					for(int i = 0; i < nbod;i++)
						acc[c][i] = 0, jerk[c][i] = 0;
				// Update Acc Jerk
				//Unroller<0,nbod*nbod>::step(bind<void>(accjerk_updater<nbod>,_1,sysref,pos,vel,acc,jerk));
				#pragma unroll
				for(int ij = nbod*nbod-1; ij >=0 ;ij--)
					accjerk_updater<nbod>(ij,sysref,pos,vel,acc,jerk);

				// Correct
				//Unroller<0,nbod>::step(bind<void>(corrector<nbod>,_1,0,sysref,pos,vel,acc,jerk,acc_old,jerk_old,h));
				for(int c = 0; c < 3 ; c++)
					#pragma unroll
					for(int i = 0; i < nbod;i++)
					corrector<nbod>(i,c,sysref,pos,vel,acc,jerk,acc_old,jerk_old,h);

				// store position and velocities in memory
				for(int c = 0; c < 3 ; c++)
					#pragma unroll
					for(int i = 0; i < nbod;i++)
						sysref[i].p(c) = pos[c][i], sysref[i].v(c) = vel[c][i];
			}

		}


		/*!
                 *  \brief Advance the system - this function must advance the system sys by one timestep, making sure that T does not exceed Tend.
		 *
		 * This function MUST return the new time of the system.
		 * This function MUST also call stop.test_body() for every body
		 * in the system, after that body has been advanced by a timestep.
		 * @tparam stop_t ...
		 * @param[in,out] ens ensemble
		 * @param[in,out] pt ...
		 * @param[in] sys system ID
		 * @param[in] T start time
		 * @param[in] Tend destination time
		 * @param[in] stop ...
		 * @param[in] stop_ts ...
		 * @param[in] step  ...
		 * @return new time of the system
		 */
		template<typename stop_t>
		__device__ double advance(ensemble &ens, thread_state_t &pt, int sys, double T, double Tend, stop_t &stop, typename stop_t::thread_state_t &stop_ts, int step)
		{
			ensemble::systemref sysref = ens[sys];
			if(T >= Tend) { return T; }
			double hh = T + this->h <= Tend ? this->h : Tend - T;

			advance_internal< 3>(sysref,hh);
			/*
			advance_internal< 4>(sysref,hh);
			advance_internal< 5>(sysref,hh);
			advance_internal< 6>(sysref,hh);
			advance_internal< 7>(sysref,hh);
			advance_internal< 8>(sysref,hh);
			advance_internal< 9>(sysref,hh);
			advance_internal<10>(sysref,hh);
*/
			return T + hh;

		}

	};

	gpu_t gpu_obj;

	/*!
         * \brief initialize temporary variables for ensemble ens. 
         *
         * This function should initialize any temporary state that is needed for integration of ens. 
	 * It will be called from gpu_generic_integrator, but only if ens.last_integrator() != this. 
         * If any temporary state exists from previous invocation of this function, it should be deallocated and the new state (re)allocated.
	 */
	void initialize(ensemble &ens)
	{
		// Here you'd initialize the object to be passed to the kernel, or
		// upload any temporary data you need to constant/texture/global memory
	}




	/*!
	 * \brief constructor 
         * 
         * Constructor will be passed the cfg object with the contents of integrator configuration file. 
         * It will be called during construction of gpu_generic_integrator. 
         * It should load any values necessary for initialization.
	 */
	prop_hermitef(const config &cfg)
	{
		if(!cfg.count("time step")) ERROR("Integrator gpu_hermitef requires a timestep ('time step' keyword in the config file).");
		gpu_obj.h = atof(cfg.at("time step").c_str());
	}
	/*!
         * \brief Cast operator for gpu_t.
         *
         * This operator must return the gpu_t object to be passed to integration kernel. 
         * It is called once per kernel invocation.
	 * @return gpu_t object to be passed to integration kernel.
	 */
	operator gpu_t()
	{
		return gpu_obj;
	}
};


/*!
 * \brief factory function to create an integrator 
 * 	  
 * This factory uses the gpu_generic_integrator class
 * with the propagator hermitef and the stopper stop_on_ejection
 * 
 * @param[in] cfg contains configuration data for gpu_generic_integrator
 */
extern "C" integrator *create_gpu_hermitef(const config &cfg)
{
	return new gpu_generic_integrator<stop_on_ejection, prop_hermitef>(cfg);
}

} // end namespace swarm

