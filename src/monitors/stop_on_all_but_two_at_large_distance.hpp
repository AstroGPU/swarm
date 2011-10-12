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

struct stop_on_all_but_two_at_large_distance_params {
	double rmax;
	stop_on_all_but_two_at_large_distance_params(const config &cfg)
	{
		double rmax = cfg.optional("rmax",std::numeric_limits<float>::max());
	}
};

/** Simple monitor that signals and logs when no more than two bodies are within a distance "rmax" of another body.
 *  \ingroup monitors
 */
template<class log_t>
class stop_on_all_but_two_at_large_distance {
	public:
	typedef stop_on_all_but_two_at_large_distance_params params;

	private:
	params _params;

	ensemble::SystemRef& _sys;
	log_t& _log;
	int _counter;

	public:

	GPUAPI bool operator () () { 
	//	if(_counter % 1000 == 0)
	//		lprintf(_log,"Hello %g\n", _sys.time() );

		_counter++;

		int num_body_near_origin = 0, id1 = -1, id2 = -2;
		for(int b = 0 ; b < _sys.nbod(); b ++ ){
			if(_sys.radius_squared(b) <= _params.rmax*_params.rmax ) // WARNING: Confusing function name
				{
				if(num_body_near_origin==0) id1 = b;
				if(num_body_near_origin==1) id2 = b;
				 num_body_near_origin++;
				}
		}
		if( num_body_near_origin > 2 ) return false;

		int num_body_far_from_all = 0;
		for(int b = 0 ; b < _sys.nbod(); b ++ ){
			if(_sys.radius_squared(b) <= _params.rmax*_params.rmax ) continue; // WARNING: Confusing function name
			bool is_far_from_every_body = true;
			for(int bb = 0 ; bb < _sys.nbod(); bb ++ ){		
			   if(b == bb) continue;
			   if(_sys.distance_squared_between(b,bb) < _params.rmax*_params.rmax )
					{ is_far_from_every_body = false;  break; }
				}
			if(is_far_from_every_body) 
				{
				num_body_far_from_all++;
				}
		}
		if(num_body_far_from_all+2>=_sys.nbod())
			{
			lprintf(_log, "No more than two bodies are within rmax of other bodies: _sys=%d, bod1=%d, bod2=%d, T=%lg r=%lg rmax=%lg.\n"
				, _sys.number(), id1, id2, _sys.time() , r, _params.rmax);
			log::system(_log, _sys);
			return true;
			}
		return false; 
	}

	GPUAPI stop_on_all_but_two_at_large_distance(const params& p,ensemble::SystemRef& s,log_t& l)
		:_params(p),_sys(s),_log(l),_counter(0){}
	
};

}
}
