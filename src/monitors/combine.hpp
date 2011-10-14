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
 *  \ingroup monitors
 */
template< template <class log_t> class Monitor1, 
	template <class log_t> class Monitor2 >
class combine_monitors {

public:

	template< class log_t >
	class combined {
		typedef Monitor1<log_t> monitor1_t;
		typedef Monitor2<log_t> monitor2_t;
		public:
		typedef combine_monitors_params<typename monitor1_t::params,typename monitor2_t::params > params;

		private:
		params _params;

		monitor1_t _monitor1;
		monitor2_t _monitor2;

		public:

		GPUAPI void operator () () { 
		  _monitor1();
		  _monitor2();
		}

		GPUAPI combined(const params& p,ensemble::SystemRef& s,log_t& l)
			:_params(p),_monitor1(p.p1,s,l),_monitor2(p.p2,s,l){}
	};
	
};

}

}
