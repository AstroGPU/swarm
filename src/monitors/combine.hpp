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
namespace monitors {

template< class Param1, class Param2 >
struct combine_monitors_params {
	Param1 p1;
	Param2 p2;
	combine_monitors_params(const config &cfg): p1(cfg), p2(cfg)
	{
	}
};

/** Template to allow developer to join two monitors
 *  Signal is true if either monitor returns true.
 *
 *  The order of the monitors is important. They are executed in order and
 *  in case of stoppers it can affect the execution. If the first one flags
 *  the system as disabled and the second one flags it as inactive then there is
 *  a problem.
 *
 *  It can be used in an expression to join more than two monitors
 *  Example:
 *  combine_monitors< combine_monitors< Monitor1, Monitor2> , Monitor3 >
 *  \ingroup monitors
 */
template< class log_t,  class Monitor1,  class Monitor2 >
struct combine {

	typedef Monitor1 monitor1_t;
	typedef Monitor2 monitor2_t;
	public:
	typedef combine_monitors_params<typename monitor1_t::params,typename monitor2_t::params > params;

	params _params;

	private:
	monitor1_t _monitor1;
	monitor2_t _monitor2;

	public:

	GPUAPI void operator () () { 
	  _monitor1();
	  _monitor2();
	}

	GPUAPI combine(const params& p,ensemble::SystemRef& s,log_t& l)
		:_params(p),_monitor1(p.p1,s,l),_monitor2(p.p2,s,l){}
	
};

}

}
