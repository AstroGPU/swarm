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

struct stop_on_crossing_orbit_params {
	stop_on_crossing_orbit_params(const config &cfg)
	{	}
};

double sqr(const double& d){
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
	params _p;
	ensemble::SystemRef& _sys;
	log_t& _log;
	
	public:

	/**
	 *  Auxiliary function to calculate (a ,e) for us
	 *  This calculation is also done in test_body from stop_on_ejection
	 *  But we would like to things simple and separate
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
		e = (fac>1.e-8) ? sqrt(fac) : 0.;
	}

    /**
	 * Function to check for crossing orbits of planet i and j.
	 * It assumes that planet i is closer to sun than planet j.
	 * WARNING: Only checks if pericenter of outer planet is less apocenter of inner planet 
	 * Doesn't account for pericenter directions
	 * Assumes planets ordered from closest to farthest
	 *
	 */
	GPUAPI bool check_for_crossing_orbits(const int& i, const int& j) {
	
	  double a_i, e_i, a_j, e_j;
	  calc_a_e(i, a_i, e_i);
	  calc_a_e(j, a_j, e_j);
	  
	  bool is_orbits_crossing = a_i * (1. + e_i)  >  a_j * ( 1. - e_j ) ;

	  if( is_orbits_crossing )
		lprintf(_log, "Crossing orbits detected: " 
				"sys=%d, T=%lg i=%d j=%d  a_i=%lg e_i=%lg a_j=%lg e_j=%lg.\n"
				, _sys.number(), _sys.time(),i,j, a_i,e_i, a_j, e_j);
		
	  return is_orbits_crossing;
	}

	GPUAPI void operator () () { 
		bool stopit = false;

		// Check for crossing orbits between every pair of planets
		// the smaller index planet always comes first
		for(int j = 0; j < _sys.nbod(); j++)
		  for(int i = 0; i < j; i++)
			  stopit = stopit || check_for_crossing_orbits(i, j);

		if(stopit) {
			log::system(_log, _sys);
			_sys.set_disabled();
		}
	}

	GPUAPI stop_on_crossing_orbit(const params& p,ensemble::SystemRef& s,log_t& l)
	    :_p(p),_sys(s),_log(l){}
	
};

} } // end namespace monitors :: swarm
