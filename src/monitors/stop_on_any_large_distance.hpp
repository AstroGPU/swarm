/*************************************************************************
 * Copyright (C) 2011 by Eric Ford and the Swarm-NG Development Team  *
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

/* Parameters for stop_on_any_large_distance monitor
 * deactivate_on_close_encounter (bool): 
 * log_on_close_encounter (bool): 
 * verbose_on_close_encounter (bool): 
 * rmax (real): minimum distance between bodies to trigger
 *
 * \ingroup monitors_param
 */ 
struct stop_on_any_large_distance_params {
	double rmax;
        bool deactivate_on, log_on, verbose_on;
	stop_on_any_large_distance_params(const config &cfg)
	{
		rmax = cfg.optional("rmax",std::numeric_limits<float>::max());
		deactivate_on = cfg.optional("deactivate_on_any_large_distance",false);
		log_on = cfg.optional("log_on_any_large_distance",false);
		verbose_on = cfg.optional("verbose_on_any_large_distance",false);
	}


};

/** Simple monitor that logs when any one body is separated from 
 *  both the origin and every other body by a distance of at least "rmax" 
 *  Optionally signals if "stop on rmax" is true
 *  This monitor may be useful in simple scattering experiments.
 *  WARNING: The code so that there is run time option of signaling to stop is untested.  
 *  TODO: Once it works, may want to port to other monitors.
 *  \ingroup monitors
 */
template<class log_t>
class stop_on_any_large_distance {
	public:
	typedef stop_on_any_large_distance_params params;

	private:
	params _params;

	ensemble::SystemRef& _sys;
	log_t& _log;

	public:
        GPUAPI bool is_deactivate_on() { return _params.deactivate_on; };
        GPUAPI bool is_log_on() { return _params.log_on; };
        GPUAPI bool is_verbose_on() { return _params.verbose_on; };
        GPUAPI bool is_any_on() { return is_deactivate_on() || is_log_on() || is_verbose_on() ; }

  //	GPUAPI void operator () () { 
	GPUAPI void operator () (int thread_in_system) 
	  if(!is_any_on()) return;

		if(thread_in_system==0)
		  {

		bool is_any_body_far_from_origin = false;
		for(int b = 0 ; b < _sys.nbod(); b ++ ){
			if(_sys.radius_squared(b) > _params.rmax * _params.rmax )
				is_any_body_far_from_origin = true;
		}
		if(!is_any_body_far_from_origin) return;
		bool need_to_log = false;
		for(int b = 0 ; b < _sys.nbod(); b ++ ){
			if(_sys.radius_squared(b) >= _params.rmax * _params.rmax ) {
				bool is_far_from_every_body = true;
				for(int bb = 0 ; bb < _sys.nbod(); bb ++ )
					if(b != bb) {
						double r2 = _sys.distance_squared_between(b,bb);
						if(r2 < _params.rmax*_params.rmax )
							is_far_from_every_body = false;
					}
				if(is_far_from_every_body) {
				  if(is_verbose_on() )
					lprintf(_log, "Distance from all bodies exceeds rmax: _sys=%d, bod=%d, T=%lg r=%lg rmax=%lg.\n", _sys.number(), b, _sys.time() , sqrt(r2), _params.rmax);
				  if(is_log_on() )
				    need_to_log = true;

				  if(is_deactivate_on() )
					a_sys.set_disabled();
				}
			}
		}
		if(need_to_log)
		   log::system(_log, _sys);		  
		  }
	}

	GPUAPI stop_on_any_large_distance(const params& p,ensemble::SystemRef& s,log_t& l)
		:_params(p),_sys(s),_log(l)
#if 0
		,_counter(0)
#endif
                {}
	
};

}
}
