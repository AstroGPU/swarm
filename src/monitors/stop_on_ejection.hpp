/*************************************************************************
 * Copyright (C) 2011 by Eric Ford and the Swarm-NG Development Team     *
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

namespace swarm {
  namespace monitors {

struct stop_on_ejection_params {
	double rmax;
        bool deactivate_on, log_on, verbose_on;
	stop_on_ejection_params(const config &cfg)
	{
		rmax = cfg.optional("rmax",std::numeric_limits<float>::max());
		deactivate_on = cfg.optional("deactivate_on_ejection",false);
		log_on = cfg.optional("log_on_ejection",false);
		verbose_on = cfg.optional("verbose_on_ejection",false);
	}
};

/** Simple monitor that signals and logs when any body (other than body 0) meets all of the following criteria:
 *	1. distance from the origin by at least rmax,  
 *	2. is moving away from the origin
 *	3. is on a nearly parabolic or hyperbolic orbit from origin (neglecting mutual interactions)
 *	Note that this stopping criteria is specifically written for planetary systems with barycenter at origin.
 *  \ingroup monitors
 */
template<class log_t>
class stop_on_ejection {
	public:
	typedef stop_on_ejection_params params;

	private:
	params _params;

	ensemble::SystemRef& _sys;
	log_t& _log;

	public:

	
        GPUAPI bool is_deactivate_on() { return _params.deactivate_on; };
        GPUAPI bool is_log_on() { return _params.log_on; };
        GPUAPI bool is_verbose_on() { return _params.verbose_on; };
        GPUAPI bool is_any_on() { return is_deactivate_on() || is_log_on() || is_verbose_on() ; }

	public:

	/** \todo This function does not look right. It compares 
	 * the radius to maximum radius, but it only stops if the
	 * shape of the orbit is parabolic or hyperbolic I thought 
	 * we should stop if the planet is too far no matter what
	 *
	 */
	GPUAPI bool test_body(const int& b) {

		double x,y,z,vx,vy,vz; _sys[b].get(x,y,z,vx,vy,vz);
		double r = sqrt(_sys[b].radius_squared());  // WARNING: Deceiving function name
		if( r < _params.rmax ) return false;
		double rdotv = x*vx+y*vy+z*vz;
		if( rdotv <= 0. ) return false;
		
		bool stopit = false;
		double speed_sq = _sys[b].speed_squared();  // WARNING: Deceiving function name
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
	
	GPUAPI void operator () () { 
	  if(!is_any_on()) return;
		bool stopit = false;

		// Check each body
		for(int b = 1; b < _sys.nbod(); b++)
			stopit = stopit || test_body(b);

		if(stopit) {
		  if(is_log_on())
			log::system(_log, _sys);
		  if(is_deactivate_on())
			_sys.set_disabled();
		}

	}

	
	GPUAPI stop_on_ejection(const params& p,ensemble::SystemRef& s,log_t& l)
		:_params(p),_sys(s),_log(l){}
	
};

}
}
