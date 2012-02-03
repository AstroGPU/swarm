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
#pragma once

#include <limits>

namespace swarm { namespace monitors {

/* Parameters for stop_on_crossing_orbit monitor
 * deactivate_on_crossing (bool): 
 * log_on_crossing (bool): 
 * verbose_on_crossing (bool): 
 *
 * \ingroup monitors_param
 *  \ingroup monitors_for_planetary_systems
 */ 
struct stop_on_crossing_orbit_params {
  bool deactivate_on, log_on, verbose_on;
  stop_on_crossing_orbit_params(const config &cfg)
	{
	  deactivate_on = cfg.optional("deactivate_on_crossing",false);
	  log_on = cfg.optional("log_on_crossing",false);
	  verbose_on = cfg.optional("verbose_on_crossing",false);
	}
};

GPUAPI double sqr(const double& d){
  return d*d;
}

/** Stopping monitor to detect crossing orbits for planets 
 *  WARNING:  This only tests for potential orbit crossing and makes assumptions about planet ordering
 *  \ingroup monitors
 */
template<class log_t>
class stop_on_crossing_orbit {
	public:
	typedef stop_on_crossing_orbit_params params;

	private:
	params _params;
	ensemble::SystemRef& _sys;
	log_t& _log;
        bool _triggered;

	public:
        GPUAPI bool is_deactivate_on() { return _params.deactivate_on; };
        GPUAPI bool is_log_on() { return _params.log_on; };
        GPUAPI bool is_verbose_on() { return _params.verbose_on; };
        GPUAPI bool is_any_on() { return is_deactivate_on() || is_log_on() || is_verbose_on() ; }

	/**
	 *  Auxiliary function to calculate (a ,e) for us
	 *  This calculation is also done in test_body from stop_on_ejection
	 *  But we would like to keep things simple and separate
	 */
	GPUAPI void calc_a_e(const int& b,double& a,double& e) {
		double x,y,z,vx,vy,vz; _sys[b].get(x,y,z,vx,vy,vz);
		// h2 = ||pos X vel||^2
		double h2 = sqr(y*vz-z*vy) + sqr(z*vx-x*vz) + sqr(x*vy-y*vx);
		double _GM = _sys[b].mass(); 
		double r = _sys[b].radius(), sp = _sys[b].speed();
		double energy = sp*0.5-_GM/r;

		a = -0.5*_GM/energy;
		double fac = 1.-h2/(_GM*a);
		e = (fac>1.e-8) ? sqrtf(fac) : 0.;
	}

    /**
	 * Function to check for crossing orbits of planet i and j.
ssssss	 * WARNING: Only checks if pericenter of outer planet is less apocenter of inner planet 
	 * Doesn't account for pericenter directions
	 * Assumes planets ordered from closest to farthest
	 *
	 */
	GPUAPI bool check_for_crossing_orbits(const int& i, const int& j) {
	
	  double a_i, e_i, a_j, e_j;
	  calc_a_e(i, a_i, e_i);
	  calc_a_e(j, a_j, e_j);
	  
	  bool is_orbits_crossing;
	  if(a_i<=a_j)
	     is_orbits_crossing = a_i * (1. + e_i)  >  a_j * ( 1. - e_j ) ;
	  else
	     is_orbits_crossing = a_i * (1. + e_i)  <  a_j * ( 1. - e_j ) ;

	  if( is_orbits_crossing && is_verbose_on())
		lprintf(_log, "Crossing orbits detected: " 
				"sys=%d, T=%lg i=%d j=%d  a_i=%lg e_i=%lg a_j=%lg e_j=%lg.\n"
				, _sys.number(), _sys.time(),i,j, a_i,e_i, a_j, e_j);
		
	  return is_orbits_crossing;
	}

#if 1 
  /// Working on standardized framework for monitors to deal with integrations in non-standard coordinate systems
	GPUAPI bool test () { 
	  _triggered = false;
		// Check for crossing orbits between every pair of planets
		for(int j = 2; j < _sys.nbod(); j++)
		  for(int i = 1; i < j; i++)
			  _triggered = _triggered || check_for_crossing_orbits(i, j);
		return _triggered;
	}
#endif

    //	GPUAPI void operator () () { 
  GPUAPI void operator () (int thread_in_system) {
	  if(!is_any_on()) return;

	  if(thread_in_system==0)
	    {
	  _triggered = test();
#if 0
	  _triggered = false;
		// Check for crossing orbits between every pair of planets
		for(int j = 2; j < _sys.nbod(); j++)
		  for(int i = 1; i < j; i++)
			  _triggered = _triggered || check_for_crossing_orbits(i, j);
#endif

		if(_triggered) {
		  if(is_log_on())
			log::system(_log, _sys);
		  if(is_deactivate_on())
			_sys.set_disabled();
		}
	    }
	}

        GPUAPI bool needs_std_coord_always () 
        {  return false; }

        GPUAPI bool needs_std_coord_now () 
        {  return false; }

	GPUAPI bool needs_to_log_system () 
        {  return (_triggered && is_log_on()); }

	GPUAPI bool needs_to_set_state () 
        {  return (_triggered && is_deactivate_on()); }


	GPUAPI stop_on_crossing_orbit(const params& p,ensemble::SystemRef& s,log_t& l)
	    :_params(p),_sys(s),_log(l){}
	
};

} } // end namespace monitors :: swarm
