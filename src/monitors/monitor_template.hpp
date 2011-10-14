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

	GPUAPI void operator () () { 
		// _sys.set_inactive();
		//_sys.set_disabled();
		// log::system(_log,_sys);
	}

	GPUAPI monitor_template(const params& p,ensemble::SystemRef& s,log_t& l)
		:_params(p),_sys(s),_log(l){}
	
};

}

}
