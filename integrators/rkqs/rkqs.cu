/*************************************************************************
 * Copyright (C) 2010 by Eric Ford    and the Swarm-NG Development Team  *
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

/*! \file rkqs.cu
 *  \brief declares prop_rkqs for use with gpu_generic_integrator
 * 
 *  based on GSL's 
 */

#include "swarm.h"
#include "../../src/swarmlog.h"
#include "rkqs.h"
#include "swarmlog.h"

/// namespace for Swarm-NG library
namespace swarm {

	// If want to declare __shared__ do it here (just like const variables for ensembles)

	/*!  
	 *  \brief propagator class for RKQS integrator on GPU: Advance the system by one time step.
	 *
	 *  CPU state and interface. Will be instantiated on construction of gpu_generic_integrator object. 
	 *  Keep any data that need to reside on the CPU here.
	 */
	struct prop_rkqs
	{
		/*! 
		 * \brief GPU state and interface (per-grid). Will be passed as an argument to integration kernel. 
		 *
		 * Any per-block read-only variables should be members of this structure.
		 */
		struct gpu_t
		{
			//! per-block variables, time step
			//! Initial time step (cfg)
			double hinit;
			//! Smallest time step allowed (cfg)
			double hmin;
			//! Largest time step allowed (cfg)
			double hmax;
			//! Error tolerance
			double error_tolerance;
		
			// Cash-Karp constants From GSL
			static const int   integrator_order = 5;
			//! Value used as power in formula to produce larger time step
			static const float step_grow_power = -1./(integrator_order+1.);
			//! Value used as power in formula to produce smaller time step
			static const float step_shrink_power = -1./integrator_order;
			//! Safety factor to prevent extreme changes in time step
			static const float step_guess_safety_factor = 0.9;
			//! Maximum growth of step size allowed at a time
			static const float step_grow_max_factor = 5.0; 
			//! Maximum shrinkage of step size allowed at a time
			static const float step_shrink_min_factor = 0.2; 

			static const int rk_order = 6;


			/*!
			 * \brief  GPU per-thread state and interface. Will be instantiated in the integration kernel. 
			 *
			 * Any per-thread variables should be members of this structure.
			 */
			struct thread_state_t
			{
				//! Current time step that it going to be used for next advance call. It will be updated by each advance call.
				double hcurrent;
				thread_state_t(const gpu_t &H, ensemble &ens, const int sys, double T, double Tend) : hcurrent(H.hinit){}
			};


				/**!
		 * helper function for acc_updater that operates on each component
		 * it gets scalar part of acceleration as input and calculates one component of
		 * acceleration at a time
         *
		 */
		template<int nbod>
		__device__ static void acc_updater_component(int i,int c
				,double dx[3],double scalar,double (&acc)[3][nbod]){
			acc[c][i] += dx[c]* scalar;
		}

		inline __device__ static double inner_product(const double a[3],const double b[3]){
			return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
		}

		inline __device__ static double sqr(const double &a){
			return a*a;
		}

		template<int nbod>
		__device__ static void acc_updater(int ij,ensemble::systemref sys
				,double (&pos)[3][nbod],double (&acc)[3][nbod]){
			const int i = ij/nbod, j = ij%nbod;
			if(i < j){

				double dx[3] =  { pos[0][j]-pos[0][i],pos[1][j]-pos[1][i],pos[2][j]-pos[2][i]};

				// computing scalar part of the acceleration
				double r2 =  dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2] ;
				double rinv = rsqrt(r2)  / r2;

				// vectorized part
				const double scalar_i = +rinv*sys[j].mass();
				acc_updater_component<nbod>(i,0,dx,scalar_i,acc);
				acc_updater_component<nbod>(i,1,dx,scalar_i,acc);
				acc_updater_component<nbod>(i,2,dx,scalar_i,acc);

				const double scalar_j = -rinv*sys[i].mass();
				acc_updater_component<nbod>(j,0,dx,scalar_j,acc);
				acc_updater_component<nbod>(j,1,dx,scalar_j,acc);
				acc_updater_component<nbod>(j,2,dx,scalar_j,acc);
			}
		}
	

		template<int nbod,int c>
			struct rk_error {
					double (&y)[3][nbod],(&k)[rk_order][3][nbod];
					const double (&b)[rk_order],h;
					rk_error(double (&y)[3][nbod],double (&k)[rk_order][3][nbod], const double (&b)[rk_order],double h)
						:y(y),k(k),b(b),h(h){}
					__device__ void operator()(int i)const{
						y[c][i] = h * (b[0]*k[0][c][i] + b[2]*k[2][c][i] + b[3]*k[3][c][i] 
								+ b[4]*k[4][c][i] + b[5]*k[5][c][i] );
					}

			};
			template<int nbod,int step, int c>
				struct rk_step {
					double (&y)[3][nbod],(&ytmp)[3][nbod],(&k)[rk_order][3][nbod];
					const double (&b)[step],h;
					rk_step(double (&y)[3][nbod],double (&ytmp)[3][nbod], double (&k)[rk_order][3][nbod], const double (&b)[step], double h): y(y),ytmp(ytmp),k(k),b(b),h(h){}
					__device__ void operator()(int i)const{
						double s = 0;
						#pragma unroll
						for(int j=0;j<step;j++)
							s += b[j] * k[j][c][i];
						ytmp[c][i] = y[c][i] + h * s;
					}
				};

			template<int nbod,int step>
			__device__ void rk_full_step(ensemble::systemref& sysref
					,double(&pos)[3][nbod],double(&in_pos)[3][nbod],double(&out_pos)[3][nbod]
					,double(&vel)[3][nbod],double(&out_vel)[3][nbod]
					,double(&vtmp)[rk_order][3][nbod],double (&atmp)[rk_order][3][nbod]
					,const double (&b)[step],double h){
				;
				for(int i = 0; i < nbod; i++) atmp[step-1][0][i] = 0;
				for(int i = 0; i < nbod; i++) atmp[step-1][1][i] = 0;
				for(int i = 0; i < nbod; i++) atmp[step-1][2][i] = 0;

				#pragma unroll
				for(int ij = nbod*nbod-1; ij >=0 ;ij--)
					acc_updater<nbod>(ij,sysref,in_pos,atmp[step-1]);
				// positions
				Unroller<0,nbod>::step(rk_step<nbod,step,0>(pos,out_pos,vtmp,b,h));
				Unroller<0,nbod>::step(rk_step<nbod,step,1>(pos,out_pos,vtmp,b,h));
				Unroller<0,nbod>::step(rk_step<nbod,step,2>(pos,out_pos,vtmp,b,h));
				// velocities
				Unroller<0,nbod>::step(rk_step<nbod,step,0>(vel,out_vel,atmp,b,h));
				Unroller<0,nbod>::step(rk_step<nbod,step,1>(vel,out_vel,atmp,b,h));
				Unroller<0,nbod>::step(rk_step<nbod,step,2>(vel,out_vel,atmp,b,h));
			}

			//! Sub stepper based on rkck (just tries a step and stores into ytmp) 
			// returns error
			template<int nbod>
			__device__ void substep(ensemble::systemref& sysref,double(&pos)[3][nbod], double(&vel)[3][nbod], double h, double& error_estimate)
			{
				if(sysref.nbod() == nbod ) {
					//		const double ah[] = { 1.0 / 5.0, 0.3, 3.0 / 5.0, 1.0, 7.0 / 8.0 };
					const double b2[] = { 1.0 / 5.0 };
					const double b3[] = { 3.0 / 40.0, 9.0 / 40.0 };
					const double b4[] = { 0.3, -0.9, 1.2 };
					const double b5[] = { -11.0 / 54.0, 2.5, -70.0 / 27.0, 35.0 / 27.0 };
					const double b6[] = { 1631.0 / 55296.0, 175.0 / 512.0, 575.0 / 13824.0, 44275.0 / 110592.0, 253.0 / 4096.0 };
					const double c1 = 37.0 / 378.0;
					const double c3 = 250.0 / 621.0;
					const double c4 = 125.0 / 594.0;
					const double c6 = 512.0 / 1771.0;

					const double cc[] = { c1, 0, c3, c4, 0, c6 };

					const double ec[] = 
						{ 
							// the first value is the same as c1, above 
							37.0 / 378.0 - 2825.0 / 27648.0,
							0.0,
							// the first value is the same as c3, above 
							250.0 / 621.0 - 18575.0 / 48384.0,
							// the first value is the same as c4, above 
							125.0 / 594.0 - 13525.0 / 55296.0,
							-277.00 / 14336.0,
							// the first value is the same as c6, above 
							512.0 / 1771.0 - 0.25 
						};



					typedef double bodarray_t[3][nbod];
					bodarray_t atmp[rk_order],ptmp,vtmp[rk_order];

					// initialize vtmp
					for(int i = 0; i < nbod; i++)
						vtmp[0][0][i] = vel[0][i],vtmp[0][1][i] = vel[1][i],vtmp[0][2][i] = vel[2][i];
					
					rk_full_step<nbod,1>(sysref,pos,pos ,ptmp,vel,vtmp[1],vtmp,atmp,b2,h);
					rk_full_step<nbod,2>(sysref,pos,ptmp,ptmp,vel,vtmp[2],vtmp,atmp,b3,h);
					rk_full_step<nbod,3>(sysref,pos,ptmp,ptmp,vel,vtmp[3],vtmp,atmp,b4,h);
					rk_full_step<nbod,4>(sysref,pos,ptmp,ptmp,vel,vtmp[4],vtmp,atmp,b5,h);
					rk_full_step<nbod,5>(sysref,pos,ptmp,ptmp,vel,vtmp[5],vtmp,atmp,b6,h);

					rk_full_step<nbod,6>(sysref,pos,ptmp,pos ,vel,vel    ,vtmp,atmp,cc,h);

					bodarray_t perr,verr;
					Unroller<0,nbod>::step(rk_error<nbod,0>(perr,vtmp,ec,h));
					Unroller<0,nbod>::step(rk_error<nbod,1>(perr,vtmp,ec,h));
					Unroller<0,nbod>::step(rk_error<nbod,2>(perr,vtmp,ec,h));
					Unroller<0,nbod>::step(rk_error<nbod,0>(verr,atmp,ec,h));
					Unroller<0,nbod>::step(rk_error<nbod,1>(verr,atmp,ec,h));
					Unroller<0,nbod>::step(rk_error<nbod,2>(verr,atmp,ec,h));

					error_estimate = max(calculate_error(pos,perr),calculate_error(vel,verr));
				}	
			}

			/// use error estimates for each coordinate to estimate the "scaled error" for the trial step
			/// based on scaling in John Chamber's mercury
			/// should we cache the scale factors (mercury does) or is recomputing better on GPU ?
			template<int nbod>
			static __device__ double calculate_error(double (&y)[3][nbod],double(&e)[3][nbod]) 
			{
				double errmax = 0.;
				for(int i=0;i<nbod;i++)
				{
					double mag = y[0][i]*y[0][i] + y[1][i]*y[1][i] + y[2][i]*y[2][i] ;
					double err = e[0][i]*e[0][i] + e[1][i]*e[1][i] + e[2][i]*e[2][i] ;

					if(err/mag>errmax) errmax = err/mag;
				}
				return errmax;
			}

			//! Utility function to calculate new time step based on error
			__device__ inline double calculate_h(double h,double error){
				double step_guess_power = (error<1.) ? step_grow_power : step_shrink_power;
				/// factor of 0.5 below due to use of squares in calculate_error, should we change to match gsl?

				/// gsl uses 1.1, but that seems dangerous, any reason we shouldn't use 1?
				double step_change_factor = ((error<0.5)||(error>1.0)) ? step_guess_safety_factor*pow(error,0.5*step_guess_power) : 1.0;

				if(error>1.)
					return max( h * max(step_change_factor,step_shrink_min_factor), hmin);
				else
					return min( h * max(min(step_change_factor,step_grow_max_factor),1.0), hmax);
			}

			/*!
			 *  \brief Advance the system - this function must advance the system sys by one timestep, making sure that T does not exceed Tend.
			 *
			 * This function MUST return the new time of the system.
			 * This function MUST also call stop.test_body() for every body
			 * in the system, after that body has been advanced by a timestep.
			 *
			 * see RKCK implementation GSL 
			 *
			 * this implementation uses RKCK as substep and updates time step based on error
			 * estimates from RKCK algorithm. Time step is stored in thread_state_t and is 
			 * updated after each time step.
			 * 
			 *
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
						if(T >= Tend) { return T; }
						// reload time step from thread state (last used time step)
						double h = T + pt.hcurrent <= Tend ? pt.hcurrent : Tend - T;
						double Ttmp = T + h;
						ensemble::systemref sysref = ens[sys];
						bool accepted = true; // the step is accepted by default to avoid infinite loop
						advance_internal<3>(sysref,h,accepted);
						advance_internal<4>(sysref,h,accepted);
						// save new timestep in thread state
						pt.hcurrent = h;
						if(accepted){
							for(int bod = 0; bod < ens.nbod(); bod++)
							{
								stop.test_body(stop_ts, ens, sys, bod, Ttmp, ens.x(sys,bod),ens.y(sys,bod),ens.z(sys,bod),ens.vx(sys,bod),ens.vy(sys,bod),ens.vz(sys,bod));
							}
							return Ttmp;
						}else
							return T;
			}

			template<int nbod>
				__device__ void advance_internal(ensemble::systemref &sysref, double & h, bool& accept_step)
				{
					if(sysref.nbod() == nbod){ 

						double pos[3][nbod],vel[3][nbod];

						#pragma unroll
						for(int i = 0; i < nbod;i++)
							pos[0][i] = sysref[i].p(0), vel[0][i] = sysref[i].v(0);
						#pragma unroll
						for(int i = 0; i < nbod;i++)
							pos[1][i] = sysref[i].p(1), vel[1][i] = sysref[i].v(1);
						#pragma unroll
						for(int i = 0; i < nbod;i++)
							pos[2][i] = sysref[i].p(2), vel[2][i] = sysref[i].v(2);

						double error_estimate = 0;
						substep<nbod>(sysref,pos,vel,h,error_estimate);
						double error = error_estimate/error_tolerance; 
						// Decision on what to is based on error, if error<1, it means error is 
						// in range and maybe we can grow time step.
						// Otherwise, calculations are erroneous and we recalculate using smaller
						// time step.
						double new_h = calculate_h(h,error);
						accept_step = (error < 1.0) || (abs(new_h - h) < 1e-10);
						// calculate new time step

						h = new_h;


						if(accept_step)
						{
							// accept the step
							#pragma unroll
							for(int i = 0; i < nbod;i++)
								sysref[i].p(0) = pos[0][i], sysref[i].v(0) = vel[0][i];
							#pragma unroll
							for(int i = 0; i < nbod;i++)
								sysref[i].p(1) = pos[1][i], sysref[i].v(1) = vel[1][i];
							#pragma unroll
							for(int i = 0; i < nbod;i++)
								sysref[i].p(2) = pos[2][i], sysref[i].v(2) = vel[2][i];
							// Store results to global memory. Invoke stopper.
						}
					}
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
		prop_rkqs(const config &cfg)
		{
			if(!cfg.count("initial time step factor")) ERROR("Integrator gpu_rkqs requires a timestep ('initial time step factor' keyword in the config file).");
			gpu_obj.hinit = atof(cfg.at("initial time step factor").c_str());
			if(!cfg.count("min time step")) ERROR("Integrator gpu_rkqs requires smallest time step allowed ('min time step' keyword in the config file).");
			gpu_obj.hmin = atof(cfg.at("min time step").c_str());
			if(!cfg.count("max time step")) ERROR("Integrator gpu_rkqs requires largest time step allowed ('min time step' keyword in the config file).");
			gpu_obj.hmax = atof(cfg.at("max time step").c_str());
			if(!cfg.count("error tolerance")) ERROR("Integrator gpu_rkqs requires an error tolerance  ('error tolerance' keyword in the config file).");
			gpu_obj.error_tolerance = atof(cfg.at("error tolerance").c_str());
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
	 * with the propagator rkqs and the stopper stop_on_ejection
	 * 
	 * @param[in] cfg contains configuration data for gpu_generic_integrator
	 */
	extern "C" integrator *create_gpu_rkqs(const config &cfg)
	{
		return new gpu_generic_integrator<stop_on_ejection, prop_rkqs>(cfg);
		//return new gpu_generic_integrator< stop_on_crossing_orbit_or_close_approach, prop_rkqs>(cfg);
	}

} // end namespace swarm
