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

struct stop_on_collision_param {
	double dmin_squared;
	stop_on_collision_param(const config &cfg)
	{
		if(!cfg.count("collision radius"))
		  dmin_squared = 0.;
		else
		  {
		    dmin_squared = atof(cfg.at("collision radius").c_str());
		    dmin_squared *= dmin_squared;
		  }
		
	}
};

/** Simple monitor to detect physical collisions.  
 *  Signals and logs if current separation between any two bodies is less than "collision radius".
 *  WARNING: Does not interpolate between steps
 *  TODO: Need to allow object specific collision radii or collision densities
 *  \ingroup monitors
 */
template<class log_t>
class stop_on_collision {
	public:
	typedef stop_on_collision_param params;

	private:
	params _p;

	ensemble::SystemRef& _sys;
	log_t& _log;

	int _counter;

	public:

	GPUAPI bool check_close_encounters(const int& i, const int& j){

		double d_squared = _sys.distance_squared_between(i,j);
		bool close_encounter = d_squared < _p.dmin_squared;

		if( close_encounter )
			lprintf(_log, "Collision detected: "
					"sys=%d, T=%f j=%d i=%d  d=%lg.\n"
				, _sys.number(), _sys.time(), j, i,sqrt(d_squared));

		return close_encounter;
	}

	GPUAPI bool operator () () { 
		bool stopit = false;

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

	GPUAPI stop_on_collision(const params& p,ensemble::SystemRef& s,log_t& l)
	    :_p(p),_sys(s),_log(l),_counter(0)){}
//		:_p(p),_sys(s),_log(l),_GM(_sys[0].mass()){}
	
};

}


