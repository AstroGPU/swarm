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

struct stop_on_close_encounter_param {
	double dmin;
	stop_on_close_encounter_param(const config &cfg)
	{
		if(!cfg.count("close approach"))
			dmin = 0.;
		else
			dmin = atof(cfg.at("close approach").c_str());
	}
};

/** Simple monitor to detect close encounters.
 *  Signals and logs if current separation between any two bodies (measured in mutual Hill radii) is less than "close approach".
 *  WARNING: Does not interpolate between steps
 *  TODO: Need to allow object specific collision radii or collision densities
 *  \ingroup monitors
 */
template<class log_t>
class stop_on_close_encounter {
	public:
	typedef stop_on_close_encounter_param params;

	private:
	params _p;

	ensemble::SystemRef& _sys;
	log_t& _log;

	int _counter;

	public:

	GPUAPI bool check_close_encounters(const int& i, const int& j){

		double d = _sys.distance_between(i,j);
		double _GM = _sys[b].mass();  // remove _ if ok to keep
		double rH = pow((_sys[i].mass()+_sys[j].mass())/(3.*_GM),1./3.);
		bool close_encounter = d < _p.dmin * rH;

		if( close_encounter )
			lprintf(_log, "Close apporach detected: "
					"sys=%d, T=%f j=%d i=%d  d=%lg.\n"
					, _sys.number(), _sys.time(), j, i,d);

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

	GPUAPI stop_on_close_encounter(const params& p,ensemble::SystemRef& s,log_t& l)
	    :_p(p),_sys(s),_log(l),_counter(0)){}
//		:_p(p),_sys(s),_log(l),_GM(_sys[0].mass()){}
	
};

}


