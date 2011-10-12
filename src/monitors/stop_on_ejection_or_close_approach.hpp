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

namespace swarm {

struct stop_on_ejection_or_close_approach_param {
	double rmax,dmin;
	stop_on_ejection_or_close_approach_param(const config &cfg)
	{
		if(!cfg.count("rmax"))
			rmax = 1000.;
		else
			rmax = atof(cfg.at("rmax").c_str());

		if(!cfg.count("close approach"))
			dmin = 0.;
		else
			dmin = atof(cfg.at("close approach").c_str());
	}
};

/** Stopping monitor to detect crossing orbits for planets or a close approach
 *  WARNING:  This has the potential to grow into a default monitor for stopping and/or logging.
 *            But I'll need to clearn this up quite a bit first.  
 *  WARNING:  This only tests for potential orbit crossing and makes assumptions about planet ordering
 *  \ingroup monitors
 */
template<class log_t>
class stop_on_ejection_or_close_approach {
	public:
	typedef stop_on_ejection_or_close_approach_param params;

	private:
	params _p;

	ensemble::SystemRef& _sys;
	log_t& _log;
	
	int _counter;
	// replaced so just use  mass of body zero
	//	double _GM;

	public:

	GPUAPI bool test_body(const int& b,double& a,double& e) {

		bool stopit = false;
		double x,y,z,vx,vy,vz; _sys[b].get(x,y,z,vx,vy,vz);
		// h2 = ||pos X vel||^2
		double h2 = sqr(y*vz-z*vy) + sqr(z*vx-x*vz) + sqr(x*vy-y*vx);
		double r = _sys[b].radius(), sp = _sys[b].speed();
		double _GM = _sys[0].mass();  // remove _ if ok to keep
		double energy = sp*0.5-_GM/r;
		double epp = energy*r/_GM;

		if( r > _p.rmax ) {
			lprintf(_log, "Distance exceeds rmax: _sys=%d, bod=%d, T=%lg r=%lg rmax=%lg.\n"
					, _sys.number(), b, _sys.time() , r, _p.rmax);
			stopit = true;
		}

		if( fabs(epp) < 1e-4 ) {
			lprintf(_log, "Orbit is parabolic: _sys=%d, bod=%d, T=%lg r=%lg energy=%lg energy*r/GM=%lg.\n"
					, _sys.number(), b, _sys.time() , r, energy, epp );
			//			stopit = true;
		}else if ( energy > 0 ){
			lprintf(_log, "Orbit is hyperbolic: _sys=%d, bod=%d, T=%lg r=%lg energy=%lg energy*r/GM=%lg.\n"
					, _sys.number(), b, _sys.time() , r, energy, epp );
			// TODO: Make sure that planet is not near another body
			// This is very unlikely to be an issue, provided that rmax
			// is set to be well beyond the initial semi-major axes
			stopit = true;
		}else {
			a = -0.5*_GM/energy;
			double fac = 1.-h2/(_GM*a);
			e = (fac>1.e-8) ? sqrt(fac) : 0.;
			stopit = false;
		}

		if(stopit) {
			lprintf(_log, "Unbound orbit detected: "
					"sys=%d, bod=%d, T=%lg, r=%lg a=%lg e=%lg.\n"
					, _sys.number(), b, _sys.time() , r, a, e);
		}

		return stopit;
	}

	GPUAPI bool check_for_crossing_orbits(
			const int& i, const double&a_i, const double & e_i
			,const int& j, const double&a_j, const double & e_j ) {

	  // WARNING: Only checks if pericenter of outer planet is less apocenter of inner planet 
	  // Doesn't account for pericenter directions
	  // Assumes planets ordered from closest to farthest
			if( a_i * (1. + e_i)  >  a_j * ( 1. - e_j ) ) {

				lprintf(_log, "Crossing orbits detected: " 
						"sys=%d, T=%lg i=%d j=%d  a_i=%lg e_i=%lg a_j=%lg e_j=%lg.\n"
						, _sys.number(), _sys.time(),i,j, a_i,e_i, a_j, e_j);
				return true;
			}else
				return false;
	}

	GPUAPI bool check_close_encounters(const int& i, const int& j){

		double d = _sys.distance_between(i,j);
		double _GM = _sys[0].mass();  // remove _ if ok to keep
		//		double rH = pow((_sys[i].mass()+_sys[j].mass())/(3.*_GM),1./3.);
		//		bool close_encounter = d < _p.dmin * rH;
		double rH3 = (_sys[i].mass()+_sys[j].mass())/(3.*_GM);
		bool close_encounter = d*d*d < _p.dmin*_p.dmin*_p.dmin * rH3;

		if( close_encounter )
			lprintf(_log, "Close apporach detected: "
					"sys=%d, T=%f j=%d i=%d  d=%lg.\n"
					, _sys.number(), _sys.time(), j, i,d);

		return close_encounter;
	}

	GPUAPI bool operator () () { 
		// WARNING: Maximum number of planet hardwired here
		double a[10],e[10];
		bool stopit = false;

		// Check each body
		for(int b = 1; b < _sys.nbod(); b++)
			stopit = stopit || test_body(b,a[b],e[b]);

		// Check for crossing orbits
		//		for(int b = 0; b < _sys.nbod()-2; b++)
		//		  stopit = stopit || check_for_crossing_orbits(b,a[b],e[b],b+1,a[b+1],e[b+1]);

		// Chcek for close encounters
		for(int b = 1; b < _sys.nbod(); b++)
			for(int d = 0; d < b; d++)
				stopit = stopit || check_close_encounters(b,d); 

		if(stopit) {
			log::system(_log, _sys);
		}

		//	if(_counter % 1000 == 0)
		//		lprintf(_log,"Hello %g\n", _sys.time() );
		_counter++;

		return stopit;
	}

	GPUAPI stop_on_ejection_or_close_approach(const params& p,ensemble::SystemRef& s,log_t& l)
	    :_p(p),_sys(s),_log(l),_counter(0)){}
//		:_p(p),_sys(s),_log(l),_GM(_sys[0].mass()){}
	
};

}


