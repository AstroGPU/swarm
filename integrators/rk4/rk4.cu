/*! \file rk4.cu
 *  \brief declares prop_rk4 for use with gpu_generic_integrator
*/

#include "swarm.h"
#include "rk4.h"
#include "swarmlog.h"

/// namespace for Swarm-NG library
namespace swarm {

/*!  
 *  \brief propagator class for RK4 integrator on GPU: Advance the system by one time step.
 *
 *  CPU state and interface. Will be instantiated on construction of gpu_generic_integrator object. 
 *  Keep any data that need to reside on the CPU here.
 */
struct prop_rk4
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
		//! per-block variables, acceleration0
		cuxDevicePtr<double, 3> aa0;
		//! per-block variables, acceleration1
		cuxDevicePtr<double, 3> aa1;
		//! per-block variables, acceleration2
		cuxDevicePtr<double, 3> aa2;

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

		/*!
                 *  \brief Advance the system - this function must advance the system sys by one timestep, making sure that T does not exceed Tend.
		 *
		 * This function MUST return the new time of the system.
		 * This function MUST also call stop.test_body() for every body
		 * in the system, after that body has been advanced by a timestep.
                 *
  		 * see RK4 implementation by http://www.artcompsci.org/kali/vol/shared_timesteps/.nbody_sh1.rb.html
		 * def rk4
                 * old_pos = pos
                 * a0 = acc(ba)
                 * pos = old_pos + vel*0.5*dt + a0*0.125*dt*dt
                 * a1 = acc(ba)
                 * pos = old_pos + vel*dt + a1*0.5*dt*dt
                 * a2 = acc(ba)
                 * pos = old_pos + vel*dt + (a0+a1*2)*(1/6.)*dt*dt
                 * vel += (a0+a1*4+a2)*(1/6.)*dt
		 * end 
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
			//double h = T + this->h <= Tend ? this->h : Tend - T;
			double h_half= h*0.5;

			/// limited to 3 bodies due to cacheing old positions and velocicties in registers (could fit 4 or 5 in 16k, more requires not cacheing v?_olds)
			double  x_old[3] = {ens.x(sys, 0), ens.x(sys, 1), ens.x(sys, 2)};
			double  y_old[3] = {ens.y(sys, 0), ens.y(sys, 1), ens.y(sys, 2)};
			double  z_old[3] = {ens.z(sys, 0), ens.z(sys, 1), ens.z(sys, 2)};
			double  vx_old[3] = {ens.vx(sys, 0), ens.vx(sys, 1), ens.vx(sys, 2)};
			double  vy_old[3] = {ens.vy(sys, 0), ens.vy(sys, 1), ens.vy(sys, 2)};
			double  vz_old[3] = {ens.vz(sys, 0), ens.vz(sys, 1), ens.vz(sys, 2)};

			// compute accelerations
			compute_acc(ens, sys, aa0);

			for(int bod = 0; bod != ens.nbod(); bod++) // starting from 1, assuming 0 is the central body
			{
				ens.x(sys,bod) = x_old[bod]+ h_half*( vx_old[bod] + aa0(sys, bod, 0) * h_half * 0.5);
				ens.y(sys,bod) = y_old[bod]+ h_half*( vy_old[bod] + aa0(sys, bod, 1) * h_half * 0.5);
				ens.z(sys,bod) = z_old[bod]+ h_half*( vz_old[bod] + aa0(sys, bod, 2) * h_half * 0.5);
			}

			// compute accelerations
			compute_acc(ens, sys, aa1);

			for(int bod = 0; bod != ens.nbod(); bod++) // starting from 1, assuming 0 is the central body
			{
				ens.x(sys,bod) = x_old[bod]+ h* (vx_old[bod] + aa1(sys, bod, 0) * h_half);
				ens.y(sys,bod) = y_old[bod]+ h* (vy_old[bod] + aa1(sys, bod, 1) * h_half);
				ens.z(sys,bod) = z_old[bod]+ h* (vz_old[bod] + aa1(sys, bod, 2) * h_half);
			}

			// compute accelerations
			compute_acc(ens, sys, aa2);

			double hby6 = h / 6.0;
			for(int bod = 0; bod != ens.nbod(); bod++) // starting from 1, assuming 0 is the central body
			{
				ens.x(sys,bod) = x_old[bod]+ h * (vx_old[bod] + (aa0(sys, bod, 0) + aa1(sys, bod, 0)*2.0) * hby6);
				ens.y(sys,bod) = y_old[bod]+ h * (vy_old[bod] + (aa0(sys, bod, 1) + aa1(sys, bod, 1)*2.0) * hby6);
				ens.z(sys,bod) = z_old[bod]+ h * (vz_old[bod] + (aa0(sys, bod, 2) + aa1(sys, bod, 2)*2.0) * hby6);
				ens.vx(sys,bod) = vx_old[bod]+ (aa0(sys, bod, 0) + 4.0* aa1(sys, bod, 0) + aa2(sys, bod, 0))* hby6;
				ens.vy(sys,bod) = vy_old[bod]+ (aa0(sys, bod, 1) + 4.0* aa1(sys, bod, 1) + aa2(sys, bod, 1))* hby6;
				ens.vz(sys,bod) = vz_old[bod]+ (aa0(sys, bod, 2) + 4.0* aa1(sys, bod, 2) + aa2(sys, bod, 2))* hby6;
				stop.test_body(stop_ts, ens, sys, bod, T+h, ens.x(sys,bod), ens.y(sys,bod), ens.z(sys,bod), ens.vx(sys,bod), ens.vy(sys,bod), ens.vz(sys,bod));
			}
		        

			return T + h;
		}

	};

	//! CPU state and interface
	cuxDeviceAutoPtr<double, 3> aa0;
	cuxDeviceAutoPtr<double, 3> aa1;
	cuxDeviceAutoPtr<double, 3> aa2;
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
		aa0.realloc(ens.nsys(), ens.nbod(), 3);
		aa1.realloc(ens.nsys(), ens.nbod(), 3);
		aa2.realloc(ens.nsys(), ens.nbod(), 3);
		gpu_obj.aa0= aa0;
		gpu_obj.aa1= aa1;
		gpu_obj.aa2= aa2;
	}

	/*!
	 * \brief constructor 
         * 
         * Constructor will be passed the cfg object with the contents of integrator configuration file. 
         * It will be called during construction of gpu_generic_integrator. 
         * It should load any values necessary for initialization.
	 */
	prop_rk4(const config &cfg)
	{
		if(!cfg.count("time step")) ERROR("Integrator gpu_rk4 requires a timestep ('time step' keyword in the config file).");
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
 * with the propagator rk4 and the stopper stop_on_ejection
 * 
 * @param[in] cfg contains configuration data for gpu_generic_integrator
 */
extern "C" integrator *create_gpu_rk4(const config &cfg)
{
	return new gpu_generic_integrator<stop_on_ejection, prop_rk4>(cfg);
}

} // end namespace swarm
