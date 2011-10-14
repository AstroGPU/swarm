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

struct stop_on_any_large_distance_or_close_encounter_params {
  double rmax,dmin ;
	stop_on_any_large_distance_or_close_encounter_params(const config &cfg)
	{
		rmax = cfg.optional("rmax",std::numeric_limits<float>::max());
		dmin = cfg.optional("close_approach", 0.);
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
class stop_on_any_large_distance_or_close_encounter {
	public:
	typedef stop_on_any_large_distance_or_close_encounter_params params;

	private:
	params _params;

	ensemble::SystemRef& _sys;
	log_t& _log;
	int _counter;

	public:

	
	public:

	GPUAPI bool test_body(const int& b) {

		double x,y,z,vx,vy,vz; _sys[b].get(x,y,z,vx,vy,vz);
		double r = sqrt(_sys[b].radius_squared());  // WARNING: Deceiving function name
		if( r < _params.rmax ) return false;
		double rdotv = x*vx+y*vy+z*vz;
		if( rdotv <= 0. ) return false;
		
		bool stopit = true;
		double speed_sq = _sys[b].speed_squared();  // WARNING: Deceiving function name
		double epp = 0.5*speed_sq*r/_sys[b].mass()-1.;
		if( fabs(epp) < 1e-4 ) {
			double energy = 0.5*speed_sq-_sys[b].mass()/r;
			lprintf(_log, "Orbit is nearly parabolic: _sys=%d, bod=%d, T=%lg r=%lg energy=%lg energy*r/GM=%lg.\n"
					, _sys.number(), b, _sys.time() , r, energy, epp );
						stopit = true;
		}else if ( epp > 0 ){
			double energy = 0.5*speed_sq-_sys[b].mass()/r;
			lprintf(_log, "Orbit is hyperbolic: _sys=%d, bod=%d, T=%lg r=%lg energy=%lg energy*r/GM=%lg.\n"
					, _sys.number(), b, _sys.time() , r, energy, epp );
						stopit = true;
		}

		return stopit;
	}
	
	GPUAPI bool check_close_encounters(const int& i, const int& j){

		double d = _sys.distance_between(i,j);
		double _GM = _sys[0].mass();  // remove _ if ok to keep
		//		double rH = pow((_sys[i].mass()+_sys[j].mass())/(3.*_GM),1./3.);
		//		bool close_encounter = d < _params.dmin * rH;
		double rH3 = (_sys[i].mass()+_sys[j].mass())/(3.*_GM);
		bool close_encounter = d*d*d < _params.dmin*_params.dmin*_params.dmin * rH3;

		if( close_encounter )
			lprintf(_log, "Close apporach detected: "
					"sys=%d, T=%f j=%d i=%d  d=%lg.\n"
					, _sys.number(), _sys.time(), j, i,d);

		return close_encounter;
	}

	GPUAPI void operator () () { 
		bool stopit = false;

		// Check each body
		for(int b = 1; b < _sys.nbod(); b++)
			stopit = stopit || test_body(b);

		// Chcek for close encounters
		for(int b = 1; b < _sys.nbod(); b++)
			for(int d = 0; d < b; d++)
				stopit = stopit || check_close_encounters(b,d); 

		if(stopit) {
			log::system(_log, _sys);
			_sys.set_disabled();
		}


	}

	
	GPUAPI stop_on_any_large_distance_or_close_encounter(const params& p,ensemble::SystemRef& s,log_t& l)
		:_params(p),_sys(s),_log(l),_counter(0){}
	
};

}
}
