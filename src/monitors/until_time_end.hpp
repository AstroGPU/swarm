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

#if 0
#include <limits>

namespace swarm {
  namespace monitors {

struct until_time_end_params {
	double time_end;
	until_time_end_params(const config &cfg)
	{
		time_end = cfg.require("time_end", 0.0);
	}
};

/** Signals true once time is equal to or greater than "time_end".
 *  Does not do any logging.
 *  Assumes integration forward in time.
 * \ingroup monitors
 *
 */
template<class log_t>
class until_time_end {
	public:
	typedef until_time_end_params params;

	private:
	params _params;

	ensemble::SystemRef& _sys;
	log_t& _log;

	public:

	GPUAPI bool operator () () { 
	  return (_sys.time() >= _params.time_end); 
	}

	GPUAPI until_time_end(const params& p,ensemble::SystemRef& s,log_t& l)
		:_params(p),_sys(s),_log(l){}
	
};

}


}
#endif
