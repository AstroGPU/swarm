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

/*! \file monitor_template.hpp
 *   \brief Defines monitor templates \ref swarm::monitors::monitor_template. 
 *
 */

#pragma once

#include <limits>

namespace swarm {

/**
 * @brief Namespace for monitors (i.e., stoppers & loggers) and their associated classes
 */
  namespace monitors {

struct monitor_template_params {
	monitor_template_params(const config &cfg)
	{
	}
};

/** Empty monitor to use as a template.  
 * Signal is always false.  Does not do any logging.
 * \ingroup monitors
 *
 */
template<class log_t>
class monitor_template {
	public:
	typedef monitor_template_params params;

	private:
	params _params;
	ensemble::SystemRef& _sys;
	log_t& _log;

	public:

        // Provide these functions, so two monitors can be combined
        GPUAPI bool is_deactivate_on() { return false; }
        GPUAPI bool is_log_on() { return false; }
        GPUAPI bool is_verbose_on() { return false; }
        GPUAPI bool is_any_on() { return is_deactivate_on() || is_log_on() || is_verbose_on() ; }
        GPUAPI bool is_condition_met () { return false; }
        GPUAPI bool need_to_log_system () 
          { return (is_log_on() && is_condition_met() ); }
        GPUAPI bool need_to_deactivate () 
          { return ( is_deactivate_on() && is_condition_met() ); }


	GPUAPI bool pass_one (int thread_in_system) 
          { return false; }
	GPUAPI int pass_two (int thread_in_system) 
          {
	    if(is_condition_met() && is_deactivate_on() && (thread_in_system==0) )
	      {  _sys.set_disabled(); }
	    return _sys.state();
	  }

        GPUAPI void operator () (const int thread_in_system) 
          { 
	    pass_one(thread_in_system);
	    pass_two(thread_in_system);
	    if(need_to_log_system() && (thread_in_system()==0) )
	      log::system(_log, _sys);
	  }

#if 0
  //	GPUAPI void operator () () { 
	GPUAPI void operator () (const int thread_in_system) { 
		// _sys.set_inactive();
		//_sys.set_disabled();
		// log::system(_log,_sys);
	}
#endif

	GPUAPI monitor_template(const params& p,ensemble::SystemRef& s,log_t& l)
		:_params(p),_sys(s),_log(l){}
	
};

}

}
