/*! \file stopper.cu
 * \brief declares stopper classes
 */

#include "swarm.h"
#include "swarmlog.h"
#include <limits>

/*

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

/// namespace for Swarm-NG library
namespace swarm {

#if 1
/// a simple stopper that triggers stopping/logging if orbits cross
/// hardwired for no more than 10 bodies
struct stop_on_crossing_orbit_or_close_approach
{
	/// GPU state and interface (per-grid) for stop_on_ejection
	struct gpu_t
	{
		float rmax;
		float dmin;


		/// GPU per-thread state and interface for stop_on_ejection
		struct thread_state_t
		{
//  WARNING: Should fix array sizes somehow
                        double a[10], e[10];
 			float GM;
			bool stop;
                        bool m_is_step_complete;

			__device__ thread_state_t(gpu_t &stop, ensemble &ens, const int sys, double T, double Tend) : stop(false), m_is_step_complete(true)
			{
                               for(int i=0;i<ens.nbod(); ++i)
                                 {
                                 a[i] = 0.;
                                 e[i] = 0.;
                                 }
                                GM = ens.mass(sys,0);
			}

                __device__ bool is_step_complete() const { return m_is_step_complete; }
                __device__ void set_step_complete() { m_is_step_complete = true; }
                __device__ void set_step_incomplete() { m_is_step_complete = false; }

		};

		/// this is called _after_ the body 'bod' has advanced a timestep.
		__device__ void test_body(thread_state_t &ts, const ensemble &ens, const int sys, const int bod, const double T, const double x, const double y, const double z, const double vx, const double vy, const double vz)
		{
			if(!ts.is_step_complete()) return;
                        if(bod==0) return; // skip sun
                        // assumes origin is barycenter, star, or jacobi center
                        double hx = y*vz-z*vy;
			double hy = z*vx-x*vz;
			double hz = x*vy-y*vx;
                        double h2 = hx*hx+hy*hy+hz*hz;
                        double hh = sqrt(h2);  
                        double r = sqrt(x*x+y*y+z*z);
                        if(r>rmax) 
{
			lprintf(dlog, "Distance exceeds rmax: sys=%d, bod=%d, T=%f r=%f rmax=%f.\n", sys, bod, T, r, rmax);
			   ts.stop = true;
       }
                 double energy = (vx*vx+vy*vy+vz*vz)*0.5-ts.GM/r;
 			if(fabs(energy*r/ts.GM)<1.e-4)
{
			lprintf(dlog, "Orbit is parabolic: sys=%d, bod=%d, T=%f r=%f energy=%f energy*r/GM=%f.\n", sys, bod, T, r, energy, energy*r/ts.GM);
		    	    ts.stop = true; //  parabola
}
			else if(energy>0) 
		{
			lprintf(dlog, "Orbit is hyperbolic: sys=%d, bod=%d, T=%f r=%f energy=%f energy*r/GM=%f.\n", sys, bod, T, r, energy, energy*r/ts.GM);
		            ts.stop = true; //  hyperbola
 }			 else
		           { // ellipse
                           ts.a[bod-1] = -0.5*ts.GM/energy;
                           double fac = 1.-h2/(ts.GM*ts.a[bod-1]);
                           ts.e[bod-1] = (fac>1e-8) ? sqrt(fac) : 0.;
			   }	
                        if(ts.stop==false) { return; }   
//			::debug_hook();
//                         ens.set_inactive(sys);
			lprintf(dlog, "Unbound orbit detected: sys=%d, bod=%d, r=%f, T=%f a=%f e=%f.\n", sys, bod, r, T,ts.a[bod-1], ts.e[bod-1]);
			log::event(dlog, swarm::log::EVT_EJECTION, T, sys, make_body_set(ens, sys, bod));
		}

		/// this is called after the entire system has completed a single timestep advance.
		__device__ bool operator ()(thread_state_t &ts, ensemble &ens, const int sys, const int step, const double T) /// should be overridden by the user
		{
			if(!ts.is_step_complete()) return false;

                        for(int i=0;i<ens.nbod()-2;++i)
                           {
                            if(ts.a[i]*(1.+ts.e[i])>ts.a[i+1]*(1.-ts.e[i+1]))
                               {
                               ts.stop =true;
                      lprintf(dlog, "Crossing orbits detected: sys=%d, T=%f i=%d i+1=%d  a_i=%f e_i=%f a_i+1=%f e_i+1=%f.\n", sys, T,i,i+1,ts.a[i], ts.e[i],ts.a[i+1],ts.e[i+1]);
                        log::event(dlog, swarm::log::EVT_EJECTION, T, sys, make_body_set(ens, sys, i));
                        log::event(dlog, swarm::log::EVT_EJECTION, T, sys, make_body_set(ens, sys, i+1));
                               }   
                            }
                        for(int i=1;i<ens.nbod();++i)
  			   {
                           float mass_i = ens.mass(sys,i);
			   for(int j=1;j<i;++j)
                              {
                              float mass_j = ens.mass(sys,j);
                              double dx = ens.x(sys,i)-ens.x(sys,j);
                              double dy = ens.y(sys,i)-ens.y(sys,j);
                              double dz = ens.z(sys,i)-ens.z(sys,j);
                              float d = sqrtf(dx*dx+dy*dy+dz*dz); 
                              double rH = pow((mass_i+mass_j)/(3.*ts.GM),1./3.);
                              if(d<dmin*rH)
                                  {
                               ts.stop =true;
                      lprintf(dlog, "Close apporach detected: sys=%d, T=%f j=%d i=%d  d=%f.\n", sys, T, j, i,d);
                        log::event(dlog, swarm::log::EVT_EJECTION, T, sys, make_body_set(ens, sys, j));
                        log::event(dlog, swarm::log::EVT_EJECTION, T, sys, make_body_set(ens, sys, i));
 				  } 
                              }
                           } 
			if(ts.stop)
			{
                               ens.set_inactive(sys);
				// store the last snapshot before going inactive
				log::system(dlog, ens, sys, T);
			}
			return ts.stop;
		}
	};

	/// CPU state and interfe for stop_on_ejection
	gpu_t gpu_obj;

	stop_on_crossing_orbit_or_close_approach(const config &cfg)
	{
		if(!cfg.count("rmax"))
		{
			gpu_obj.rmax = 0.;
		}
		else
		{
			gpu_obj.rmax = atof(cfg.at("rmax").c_str());
		}
		if(!cfg.count("close approach"))
		{
			gpu_obj.dmin = 0.;
		}
		else
		{
			gpu_obj.dmin = atof(cfg.at("close approach").c_str());
		}
	}

 	// currently empty for stop_on_ejection
	void initialize(const ensemble &ens)
	{
		// Here you'd upload any temporary data you need to constant/texture/global memory
	}

 	// return gpu data for stop_on_ejection
	operator gpu_t()
	{
		return gpu_obj;
	}
};
#endif



/// a simple stopper that checks if the distance from origin exceeds a threshold (for demonstration purposes only)
struct stop_on_ejection
{
	/// GPU state and interface (per-grid) for stop_on_ejection
	struct gpu_t
	{
		float rmax;

		/// GPU per-thread state and interface for stop_on_ejection
		struct thread_state_t
		{
			bool eject;
			bool m_is_step_complete;

			__device__ thread_state_t(gpu_t &stop, ensemble &ens, const int sys, double T, double Tend) :
				eject(false), m_is_step_complete(true)
			{
			}

                __device__ bool is_step_complete() const { return m_is_step_complete; }
                __device__ void set_step_complete() { m_is_step_complete = true; }
                __device__ void set_step_incomplete() { m_is_step_complete = false; }

		};

		/// this is called _after_ the body 'bod' has advanced a timestep.
		__device__ void test_body(thread_state_t &ts, const ensemble &ens, const int sys, const int bod, const double T, const double x, const double y, const double z, const double vx, const double vy, const double vz)
		{
			if(!ts.is_step_complete()) return;
			float r = sqrtf(x*x + y*y + z*z);
			if(r < rmax) { return; }
//			::debug_hook();
			ts.eject = true;
//			lprintf(dlog, "Ejection detected: sys=%d, bod=%d, r=%f, T=%f.\n", sys, bod, r, T);
			log::event(dlog, swarm::log::EVT_EJECTION, T, sys, make_body_set(ens, sys, bod));
		}

		/// this is called after the entire system has completed a single timestep advance.
		__device__ bool operator ()(thread_state_t &ts, ensemble &ens, const int sys, const int step, const double T) /// should be overridden by the user
		{
			if(!ts.is_step_complete()) return false;
			if(ts.eject)
			{
				// store the last snapshot before going inactive
				log::system(dlog, ens, sys, T);
			}
			return ts.eject;
		}
	};

	/// CPU state and interface for stop_on_ejection
	gpu_t gpu_obj;

	stop_on_ejection(const config &cfg)
	{
		if(!cfg.count("rmax"))
		{
			gpu_obj.rmax = std::numeric_limits<float>::max();
		}
		else
		{
			gpu_obj.rmax = atof(cfg.at("rmax").c_str());
		}
	}

 	// currently empty for stop_on_ejection
	void initialize(const ensemble &ens)
	{
		// Here you'd upload any temporary data you need to constant/texture/global memory
	}

 	// return gpu data for stop_on_ejection
	operator gpu_t()
	{
		return gpu_obj;
	}
};

} // end namespace swarm
