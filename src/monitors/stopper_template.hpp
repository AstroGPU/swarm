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


struct stopper_template_params {
	stopper_template_params(const config &cfg)
	{
	}
};

/** Empty monitor to use as a template
 * \ingroup monitors
 *
 */
template<class log_t>
class stopper_template {
	public:
	typedef stopper_template_params params;

	private:
	params _params;

	ensemble::SystemRef& _sys;
	log_t& _log;

	public:

	GPUAPI bool operator () () { 
		return false; 
	}

	GPUAPI stopper_template(const params& p,ensemble::SystemRef& s,log_t& l)
		:_params(p),_sys(s),_log(l){}
	
};

}

