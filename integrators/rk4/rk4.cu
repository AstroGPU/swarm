#include "swarm.h"
#include "rk4.h"
#include "swarmlog.h"

/*!

	\brief gpu_generic_integrator - Versatile GPU integrator template

	gpu_generic_integrator is a template class designed to make it easy to implement
	powerful GPU-based integrators by supplying a 'propagator' class (that provides
	a function to advance the ensemble by one timestep), and a 'stopper' class (that
	tests whether the integration should stop for a particular system in the
	ensemble). Given the two classes, gpu_generic_integrator acts as a driver
	routine repeatedly calling the GPU kernel until all systems in the ensemble have
	finished the integration. It does this in an optimized manner, taking care to
	compactify() the ensemble when needed to efficiently utilize GPU resources
	(NOTE: compactification not yet implemented).

	For a cannonical example of how to build an integrator using
	gpu_generic_integrator, look at the gpu_euler integrator.

	Integration loop outline (I == gpu_generic_integrator object, H == propagator
	object, stop == stopper object):

	I.integrate(ens):
		if(ens.last_integrator() != this):
			H.initialize(ens);
			stop.initialize(ens);
		do:
			gpu_integrate<<<>>>(max_step, ens, (gpu_t)H, (gpu_t)stop)	[m'threaded execution on the GPU]
				while step < max_step:
					if(T < Tend || stop()):
						ens(sys).flags |= INACTIVE
						break;
					H.advance():			[ implementation supplied by the developer ]
						foreach(bod in sys):
							advance bod
							call stop.test_body
						return new time T
			if(beneficial):
				ens.compactify();
		while(active_systems > 0)

	To build an integrator using gpu_generic_integrator, the developer must supply a
	propagator and a stopper class that conform to the following interfaces:

	// propagator class: advance the system by one time step
	//
	// CPU state and interface. Will be instantiated on construction
	// of gpu_generic_integrator object. Keep any data that needs
	// to reside on the CPU here.
	struct propagator
	{
		// GPU state and interface (per-grid). Will be passed 
		// as an argument to integration kernel. Any per-block read-only
		// variables should be members of this structure.
		struct gpu_t
		{
			// GPU per-thread state and interface. Will be instantiated
			// in the integration kernel. Any per-thread
			// variables should be members of this structure.
			struct thread_state_t
			{
				__device__ thread_state_t(const gpu_t &H, ensemble &ens, const int sys, double T, double Tend);
			};

			// Advance the system - this function must advance the system
			// sys by one timestep, making sure that T does not exceed Tend.
			// Must return the new time of the system.
			//
			// This function MUST also call stop.test_body() for every body
			// in the system, after that body has been advanced by a timestep.
			template<typename stop_t>
			__device__ double advance(ensemble &ens, thread_state_t &pt, int sys, double T, double Tend, stop_t &stop, typename stop_t::thread_state_t &stop_ts, int step)
		};

		// Constructor will be passed the cfg object with the contents of
		// integrator configuration file. It will be called during construction
		// of gpu_generic_integrator. It should load any values necessary
		// for initialization.
		propagator(const config &cfg);

		// Initialize temporary variables for ensemble ens. This function
		// should initialize any temporary state that is needed for integration
		// of ens. It will be called from gpu_generic_integrator, but only
		// if ens.last_integrator() != this. If any temporary state exists from
		// previous invocation of this function, it should be deallocated and
		// the new state (re)allocated.
		void initialize(ensemble &ens);

		// Cast operator for gpu_t. This operator must return the gpu_t object
		// to be passed to integration kernel. It is called once per kernel
		// invocation.
		operator gpu_t();
	};


	// stopper class: mark a system inactive if conditions are met
	//
	// CPU state and interface. Will be instantiated on construction
	// of gpu_generic_integrator object. Keep any data that needs
	// to reside on the CPU here.
	struct stopper
	{
		// GPU state and interface (per-grid). Will be passed 
		// as an argument to integration kernel. Any per-block read-only
		// variables should be members of this structure.
		struct gpu_t
		{
			// GPU per-thread state and interface. Will be instantiated
			// in the integration kernel. Any per-thread
			// variables should be members of this structure.
			struct thread_state_t
			{
				__device__ thread_state_t(gpu_t &stop, ensemble &ens, const int sys, double T, double Tend);
			};

			// test any per-body stopping criteria for body (sys,bod). If 
			// your stopping criterion only depends on (x,v), test for it 
			// here. This will save you the unnecessary memory accesses 
			// that would otherwise be made if the test was made from 
			// operator().
			//
			// Called _after_ the body 'bod' has advanced a timestep.
			//
			// Note: you must internally store the result of your test,
			// and use/return it in subsequent call to operator().
			//
			__device__ void test_body(thread_state_t &ts, ensemble &ens, int sys, int bod, double T, double x, double y, double z, double vx, double vy, double vz);

			// Called after a system sys has been advanced by a timestep.
			// Must return true if the system sys is to be flagged as
			// INACTIVE (thus stopping further integration)
			__device__ bool operator ()(thread_state_t &ts, ensemble &ens, int sys, int step, double T);
		};

		// Constructor will be passed the cfg object with the contents of
		// integrator configuration file. It will be called during construction
		// of gpu_generic_integrator. It should load any values necessary
		// for initialization.
		stopper(const config &cfg);

		// Initialize temporary variables for ensemble ens. This function
		// should initialize any temporary state that is needed for integration
		// of ens. It will be called from gpu_generic_integrator, but only
		// if ens.last_integrator() != this. If any temporary state exists from
		// previous invocation of this function, it should be deallocated and
		// the new state (re)allocated.
		void initialize(ensemble &ens);

		// Cast operator for gpu_t. This operator must return the gpu_t object
		// to be passed to integration kernel. It is called once per kernel
		// invocation.
		operator gpu_t();
	};

*/

namespace swarm {


/*!  
 *  \brief propagator class for RK4 integrator on GPU: Advance the system by one time step.
 *
 *  CPU state and interface. Will be instantiated on construction of gpu_generic_integrator object. 
 *  Keep any data that needs to reside on the CPU here.
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
                 * vel += (a0_a1*4+a2)*(1/6.)*dt
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
				ens.x(sys,bod) = x_old[bod]+ vx_old[bod] * h_half + aa0(sys, bod, 0) * h_half * h_half * 0.5;
				ens.y(sys,bod) = y_old[bod]+ vy_old[bod] * h_half + aa0(sys, bod, 1) * h_half * h_half * 0.5;
				ens.z(sys,bod) = z_old[bod]+ vz_old[bod] * h_half + aa0(sys, bod, 2) * h_half * h_half * 0.5;
			}

			// compute accelerations
			compute_acc(ens, sys, aa1);

			for(int bod = 0; bod != ens.nbod(); bod++) // starting from 1, assuming 0 is the central body
			{
				ens.x(sys,bod) = x_old[bod]+ vx_old[bod] * h + aa1(sys, bod, 0) * h_half * h;
				ens.y(sys,bod) = y_old[bod]+ vy_old[bod] * h + aa1(sys, bod, 1) * h_half * h;
				ens.z(sys,bod) = z_old[bod]+ vz_old[bod] * h + aa1(sys, bod, 2) * h_half * h;
			}

			// compute accelerations
			compute_acc(ens, sys, aa2);

			for(int bod = 0; bod != ens.nbod(); bod++) // starting from 1, assuming 0 is the central body
			{
				ens.x(sys,bod) = x_old[bod]+ vx_old[bod] * h + (aa0(sys, bod, 0) + aa1(sys, bod, 0)*2.0) * h * h / 6.0;
				ens.y(sys,bod) = y_old[bod]+ vy_old[bod] * h + (aa0(sys, bod, 1) + aa1(sys, bod, 1)*2.0) * h * h / 6.0;
				ens.z(sys,bod) = z_old[bod]+ vz_old[bod] * h + (aa0(sys, bod, 2) + aa1(sys, bod, 2)*2.0) * h * h / 6.0;
				ens.vx(sys,bod) = vx_old[bod]+ (aa0(sys, bod, 0) + 4.0* aa1(sys, bod, 0) + aa2(sys, bod, 0))* h/6.0;
				ens.vy(sys,bod) = vy_old[bod]+ (aa0(sys, bod, 1) + 4.0* aa1(sys, bod, 1) + aa2(sys, bod, 1))* h/6.0;
				ens.vz(sys,bod) = vz_old[bod]+ (aa0(sys, bod, 2) + 4.0* aa1(sys, bod, 2) + aa2(sys, bod, 2))* h/6.0;
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
		if(!cfg.count("h")) ERROR("Integrator gpu_rk4 needs a timestep ('h' keyword in the config file).");
		gpu_obj.h = atof(cfg.at("h").c_str());
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


// factory
extern "C" integrator *create_gpu_rk4(const config &cfg)
{
	return new gpu_generic_integrator<stop_on_ejection, prop_rk4>(cfg);
}

} // end namespace swarm
