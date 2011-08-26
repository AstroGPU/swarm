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


namespace swarm {


struct log_time_interval_params {
	double time_interval;
	log_time_interval_params(const config &cfg)
	{
		time_interval = cfg.require("log interval", 0.0);
	}
};

template<class log_t>
class log_time_interval {
	public:
	typedef log_time_interval_params params;

	private:
	params _params;

	ensemble::SystemRef& _sys;
	double _next_log_time;
	log_t& _log;

	public:

	GPUAPI bool operator () () { 
		if( _sys.time() > _next_log_time ) {
			log::system(_log, _sys );
			_next_log_time += _params.time_interval; 
			//lprintf(_log,"Logging at %lg\n",_sys.time());
		}
		return false;
	}

	GPUAPI log_time_interval(const params& p,ensemble::SystemRef& s,log_t& l)
		:_params(p),_sys(s),_log(l),_next_log_time(s.time()){}
	
};

}


