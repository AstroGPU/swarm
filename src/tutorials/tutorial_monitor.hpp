/*************************************************************************
 * Copyright (C) 2009-2010 by Eric Ford & the Swarm-NG Development Team  *
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

/** \file tutorial_monitor.hpp
 *  \brief A tutorial for making a monitor.
 *
 */

/*
 *  This is a simple tutorial used in doxygen pages
 *  should go through program2doxygen before it
 *  can be used by doxygen.
 *  
 *
 */
// \page TutorialMonitor Tutorial for making a monitor (stopper/logger)
//
#pragma once

#include <limits>

namespace swarm { namespace monitors {

// Parameters for TutorialMonitor monitor
// deactivate_on_ejection (bool): 
// log_on_ejection (bool): 
// verbose_on_ejection (bool): 
// rmax (real): minimum distance to check for ejections
//
struct TutorialMonitor_params {
	double rmax;
        bool deactivate_on, log_on, verbose_on;
	TutorialMonitor_params(const config &cfg)
	{
		rmax = cfg.optional("rmax",std::numeric_limits<float>::max());
		deactivate_on = cfg.optional("deactivate_on_ejection",false);
		log_on = cfg.optional("log_on_ejection",false);
		verbose_on = cfg.optional("verbose_on_ejection",false);
	}
};

// Simple monitor that signals and logs when any body (other than body 0) meets all of the following criteria:
//	1. distance from the origin by at least rmax,  
//	2. is moving away from the origin
//	3. is on a nearly parabolic or hyperbolic orbit from origin (neglecting mutual interactions)
//	Note that this stopping criteria is specifically written for planetary systems with barycenter or star at origin.
template<class log_t>
class TutorialMonitor {
	public:
	typedef TutorialMonitor_params params;

	private:
	params _params;
        bool condition_met;

	ensemble::SystemRef& _sys;
	log_t& _log;

	public:

        GPUAPI bool is_deactivate_on() { return _params.deactivate_on; };
        GPUAPI bool is_log_on() { return _params.log_on; };
        GPUAPI bool is_verbose_on() { return _params.verbose_on; };
        GPUAPI bool is_any_on() { return is_deactivate_on() || is_log_on() || is_verbose_on() ; }
        GPUAPI bool is_condition_met () { return ( condition_met ); }
        GPUAPI bool need_to_log_system () 
          { return (is_log_on() && is_condition_met() ); }
        GPUAPI bool need_to_deactivate () 
          { return ( is_deactivate_on() && is_condition_met() ); }

        GPUAPI void log_system()  {  log::system(_log, _sys);  }

		template<class T>
		static GENERIC int thread_per_system(T compile_time_param){
			return 1;
		}

		template<class T>
		static GENERIC int shmem_per_system(T compile_time_param) {
			 return 0;
		}

	public:

	// First, we test whether body b is far from the origin 
	// (which may be COM or star) and moving farther away.  
	// If both of those are passed, then we check that the
	// orbit is parabolic or hyperbolic.
	//
	GPUAPI bool test_body(const int& b) {

		double x,y,z,vx,vy,vz; _sys[b].get(x,y,z,vx,vy,vz);
		double r = _sys[b].radius();  
		if( r < _params.rmax ) return false;
		double rdotv = x*vx+y*vy+z*vz;
		if( rdotv <= 0. ) return false;
		
		bool stopit = false;
		double speed_sq = _sys[b].speed_squared();  
		double epp = 0.5*speed_sq*r/_sys[b].mass()-1.;
		if( fabs(epp) < 1e-4 ) {
			double energy = 0.5*speed_sq-_sys[b].mass()/r;
			if(is_verbose_on() )
			  lprintf(_log, "Orbit is nearly parabolic: _sys=%d, bod=%d, T=%lg r=%lg energy=%lg energy*r/GM=%lg.\n"
					, _sys.number(), b, _sys.time() , r, energy, epp );
			stopit = true;
		}else if ( epp > 0 ){
			double energy = 0.5*speed_sq-_sys[b].mass()/r;
			if(is_verbose_on() )
			  lprintf(_log, "Orbit is hyperbolic: _sys=%d, bod=%d, T=%lg r=%lg energy=%lg energy*r/GM=%lg.\n"
					, _sys.number(), b, _sys.time() , r, energy, epp );
			// TODO: Make sure that planet is not near another body
			// This is very unlikely to be an issue, provided that rmax
			// is set to be well beyond the initial semi-major axes
			stopit = true;
		}

		return stopit;
	}
	
        GPUAPI void operator () (const int thread_in_system) 
          { 
	    pass_one(thread_in_system);
	    pass_two(thread_in_system);
	    if(need_to_log_system() && (thread_in_system==0) )
	      log_system();
	  }


	
	GPUAPI TutorialMonitor(const params& p,ensemble::SystemRef& s,log_t& l)
		:_params(p),_sys(s),_log(l){}

	GPUAPI bool pass_one (int thread_in_system) 
          {
	    bool need_full_test = false; 
	    condition_met = false;
	    if(is_any_on()&&(thread_in_system==0))
	      {
		// Check each body other than central star
		for(int b = 1; b < _sys.nbod(); b++)
		  condition_met = condition_met || test_body(b);
		
		if( condition_met && is_log_on() )
		  {  need_full_test = true;  }


	      }
	    return need_full_test;
	  }
	    

	GPUAPI int pass_two (int thread_in_system) 
          {
	    if(is_condition_met() && is_deactivate_on() &&(thread_in_system==0) )
	      {  _sys.set_disabled(); }
	    return _sys.state();
	  }

	
};

} } // close namespace monitor :: swarm
