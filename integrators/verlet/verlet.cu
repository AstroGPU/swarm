#include "swarm.h"
#include "verlet.h"
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

	For a canonical example of how to build an integrator using
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
	// of gpu_generic_integrator object. Keep any data that need
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
	// of gpu_generic_integrator object. Keep any data that need
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

/// namespace for Swarm-NG library
namespace swarm {


/*!  
 *  \brief propagator class for verlet integrator on GPU: Advance the system by one time step.
 *
 *  CPU state and interface. Will be instantiated on construction of gpu_generic_integrator object. 
 *  Keep any data that need to reside on the CPU here.
 */
struct prop_verlet
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
		//! per-block variables, acceleration
		cuxDevicePtr<double, 3> aa;

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

			//Step(pos);
			//CalcDerivForDrift(); -> returns velocity
			for(int bod = 0; bod != ens.nbod(); bod++) // starting from 1, assuming 0 is the central body
			{
				double  x = ens.x(sys, bod),   y = ens.y(sys, bod),   z = ens.z(sys, bod);
				double vx = ens.vx(sys, bod), vy = ens.vy(sys, bod), vz = ens.vz(sys, bod);
				x += vx * h_half;
				y += vy * h_half;
				z += vz * h_half;
				ens.x(sys, bod)  = x;   ens.y(sys, bod) = y;   ens.z(sys, bod) = z;
			}

			//CalcDerivForKick();
			// compute accelerations
			compute_acc(ens, sys, aa);

			//Step(vel)
			for(int bod = 0; bod != ens.nbod(); bod++) 
			{
				// load
				double  x = ens.x(sys, bod),   y = ens.y(sys, bod),   z = ens.z(sys, bod);
				double vx = ens.vx(sys, bod), vy = ens.vy(sys, bod), vz = ens.vz(sys, bod);
				// advance
				//x += vx * h_half;
				//y += vy * h_half;
				//z += vz * h_half;
				vx += aa(sys, bod, 0) * h_half;
				vy += aa(sys, bod, 1) * h_half;
				vz += aa(sys, bod, 2) * h_half;

				//stop.test_body(stop_ts, ens, sys, bod, T+h_half, x, y, z, vx, vy, vz);

				// store
				//ens.x(sys, bod)  = x;   ens.y(sys, bod) = y;   ens.z(sys, bod) = z;
				ens.vx(sys, bod) = vx; ens.vy(sys, bod) = vy; ens.vz(sys, bod) = vz;
			}
		        
			T = T + h_half;
			//Calculate New HalfDeltaT; HalfDeltaT =1./(2./(CalcTimeScaleFactor(mass, pos, nBodies)*DeltaTau)-1./HalfDeltaTLast);
			//CalcTimeScaleFactor(mass, pos, nBodies)
			double rinv3 = 0.;
			for ( unsigned int i=0;i<ens.nbod();++i )
			{
				//V3 xi( ens.x ( sys,i ), ens.y ( sys,i ), ens.z ( sys,i ) );
				double xi0= ens.x ( sys,i );
				double xi1= ens.y ( sys,i );
				double xi2= ens.z ( sys,i );
				
				double mass_sum = ens.mass(sys,i);
				for (int j=i+1; j < ens.nbod(); ++j) {
					//double msum = g_mass[i]+g_mass[j];
					mass_sum += ens.mass(sys,j); 
					//V3 dx(ens.x(sys,j), ens.y(sys,j), ens.z(sys,j));  dx -= xi;
					//double r2 = dx.MagnitudeSquared();
					double dx0=ens.x(sys,j) - xi0;
					double dx1=ens.y(sys,j) - xi1;
					double dx2=ens.z(sys,j) - xi2;
					double r2 = dx0*dx0 + dx1*dx1 + dx2*dx2;
					double rinv = 1./sqrt ( r2 );
					rinv *= mass_sum;
					rinv3 += rinv/r2;
				}
			}

			double h_new= 1./(1.+sqrt(rinv3)); 
			// Technically, missing a factor of 2.*M_PI, but this choice is arbitrary 
			h_new = 1./(2./(h_new*h)-1./h_half);
			h_half = h_new;

			//Step(vel)
			for(int bod = 0; bod != ens.nbod(); bod++) 
			{
				// load
				double  x = ens.x(sys, bod),   y = ens.y(sys, bod),   z = ens.z(sys, bod);
				double vx = ens.vx(sys, bod), vy = ens.vy(sys, bod), vz = ens.vz(sys, bod);
				vx += aa(sys, bod, 0) * h_half;
				vy += aa(sys, bod, 1) * h_half;
				vz += aa(sys, bod, 2) * h_half;
				x += vx * h_half;
				y += vy * h_half;
				z += vz * h_half;

				//stop.test_body(stop_ts, ens, sys, bod, T+h_half, x, y, z, vx, vy, vz);

				// store
				ens.x(sys, bod)  = x;   ens.y(sys, bod) = y;   ens.z(sys, bod) = z;
				ens.vx(sys, bod) = vx; ens.vy(sys, bod) = vy; ens.vz(sys, bod) = vz;
			}


			return T + h_half;
		}
	};

	//! CPU state and interface
	cuxDeviceAutoPtr<double, 3> aa;
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
		aa.realloc(ens.nsys(), ens.nbod(), 3);
		gpu_obj.aa = aa;
	}

	/*!
	 * \brief constructor 
         * 
         * Constructor will be passed the cfg object with the contents of integrator configuration file. 
         * It will be called during construction of gpu_generic_integrator. 
         * It should load any values necessary for initialization.
	 */
	prop_verlet(const config &cfg)
	{
		if(!cfg.count("time step")) ERROR("Integrator gpu_verlet needs a timestep ('time step' keyword in the config file).");
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


// factory
extern "C" integrator *create_gpu_verlet(const config &cfg)
{
	return new gpu_generic_integrator<stop_on_ejection, prop_verlet>(cfg);
}

} // end namespace swarm
