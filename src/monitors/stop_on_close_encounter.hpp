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

namespace swarm { namespace monitors {

struct stop_on_close_encounter_param {
	double dmin;
  bool deactivate_on, log_on, verbose_on;
	stop_on_close_encounter_param(const config &cfg)
	{
		dmin = cfg.optional("close_approach",false);
		deactivate_on = cfg.optional("deactivate_on_close_encounter",false);
		log_on = cfg.optional("log_on_close_encounter",false);
		verbose_on = cfg.optional("verbose_on_close_encounter",false);
	}
};

/** Simple monitor to detect close encounters.
 *  Signals and logs if current separation between any two bodies (measured in mutual Hill radii) is less than "close_approach".
 *  WARNING: Does not interpolate between steps
 *  TODO: Need to allow object specific collision radii or collision densities
 *  \ingroup monitors
 */
template<class log_t>
class stop_on_close_encounter {
	public:
	typedef stop_on_close_encounter_param params;

	private:
	params _params;
	ensemble::SystemRef& _sys;
	log_t& _log;


	public:

        GPUAPI bool is_deactivate_on() { return _params.deactivate_on; };
        GPUAPI bool is_log_on() { return _params.log_on; };
        GPUAPI bool is_verbose_on() { return _params.verbose_on; };
        GPUAPI bool is_any_on() { return is_deactivate_on() || is_log_on() || is_verbose_on() ; }


	GPUAPI bool check_close_encounters(const int& i, const int& j){

		double d = _sys.distance_between(i,j);
		double _GM = _sys[0].mass();  // remove _ if ok to keep
		//		double rH = pow((_sys[i].mass()+_sys[j].mass())/(3.*_GM),1./3.);
		//		bool close_encounter = d < _p.dmin * rH;
		double rH3 = (_sys[i].mass()+_sys[j].mass())/(3.*_GM);
		bool close_encounter = d*d*d < _params.dmin*_params.dmin*_params.dmin * rH3;

		if( close_encounter )
		  if(is_verbose_on() )
			lprintf(_log, "Close apporach detected: "
					"sys=%d, T=%f j=%d i=%d  d=%lg.\n"
					, _sys.number(), _sys.time(), j, i,d);

		return close_encounter;
	}

	GPUAPI void operator () () { 
	  if(!is_any_on()) return;
		bool stopit = false;

		// Chcek for close encounters
		for(int b = 1; b < _sys.nbod(); b++)
			for(int d = 0; d < b; d++)
				stopit = stopit || check_close_encounters(b,d); 

		if(stopit) {
		  if(is_log_on())
			log::system(_log, _sys);
		  if(is_deactivate_on())
			_sys.set_disabled();
		}
	}

	GPUAPI stop_on_close_encounter(const params& p,ensemble::SystemRef& s,log_t& l)
	    :_params(p),_sys(s),_log(l){}
	
};

} } // end namespace monitors :: swarm
