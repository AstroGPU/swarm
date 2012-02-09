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
  namespace monitors {

/* Parameters for log_time_interval monitor
 * log_on_interval (bool): 
 * log_interval (real): time between sucessive logging
 * \ingroup monitors_param
 */ 
struct log_time_interval_params {
	double time_interval;
        bool log_on;
	log_time_interval_params(const config &cfg)
	{
		time_interval = cfg.require("log_interval", 0.0);
		log_on = cfg.optional("log_on_interval",true);
	}
};

/** Monitor that logs the entire state of systems at periodic intervals of approximately "log_interval"
 *  Systems may be integrated for more than log interval before another log entry is written.
 *  Assumes integration results in increasing time.
 * 
 *  \ingroup monitors
 *
 *   the period of logging should be specified as "log_interval" in config file.
 */
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
        GPUAPI bool is_deactivate_on() { return false; };
        GPUAPI bool is_log_on() { return _params.log_on; };
        GPUAPI bool is_verbose_on() { return false; };
        GPUAPI bool is_any_on() { return is_deactivate_on() || is_log_on() || is_verbose_on() ; }

	GPUAPI void operator () (const int thread_in_system) { 
	  if( (thread_in_system==0) && is_log_on() )
	    {
		if(_sys.time() >= _next_log_time )  {
			log::system(_log, _sys );
			_next_log_time += _params.time_interval; 
			//lprintf(_log,"Logging at %lg\n",_sys.time());
		}
	    }
	}

	GPUAPI log_time_interval(const params& p,ensemble::SystemRef& s,log_t& l)
		:_params(p),_sys(s),_log(l),_next_log_time(s.time()){}
	
};

}


}
