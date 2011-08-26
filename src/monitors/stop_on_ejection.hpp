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


struct stop_on_ejection_params {
	double rmax_squared;
	stop_on_ejection_params(const config &cfg)
	{
		double rmax = cfg.optional("rmax",std::numeric_limits<float>::max());
		rmax_squared = rmax * rmax;
	}
};

/** Simple monitor that stops the system when a planet ejects from the system
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
	int _counter;

	public:

	GPUAPI bool operator () () { 
		for(int b = 1 ; b < _sys.nbod(); b ++ ){
			if(_sys.distance_squared_between(b,0) > _params.rmax_squared )
				return true;
		}
	//	if(_counter % 1000 == 0)
	//		lprintf(_log,"Hello %g\n", _sys.time() );

		_counter++;

		return false; 
	}

	GPUAPI stop_on_ejection(const params& p,ensemble::SystemRef& s,log_t& l)
		:_params(p),_sys(s),_log(l),_counter(0){}
	
};

}
