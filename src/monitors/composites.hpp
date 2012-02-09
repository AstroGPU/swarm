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

#include "combine.hpp"
#include "stop_on_ejection.hpp"
#include "stop_on_crossing_orbit.hpp"
#include "stop_on_close_encounter.hpp"

namespace swarm { namespace monitors {


/** Combination of stop_on_ejcetion, stop_on_close_encounter and
 * stop_on_crossing_orbit
 *
 */
template <class L> 
struct stop_on_ejection_or_close_encounter_or_crossing_orbit {
	struct params {
		typename stop_on_ejection<L>        ::params ej;
		typename stop_on_close_encounter<L> ::params ce;
		typename stop_on_crossing_orbit<L>  ::params co;
		
		params(const config& cfg)
			:ej(cfg), ce(cfg), co(cfg) {}
	};
	
	GPUAPI stop_on_ejection_or_close_encounter_or_crossing_orbit
		(const params& p,ensemble::SystemRef& s,L& l)
		: ej(p.ej,s,l), ce(p.ce,s,l), co(p.co,s,l) 	{}
	
  //	GPUAPI void operator () ()
  //         {	  ej(); ce(); co(); }

	GPUAPI void operator () (const int thread_in_system) 
         {  ej(thread_in_system); ce(thread_in_system); co(thread_in_system);  }

	
private:
	stop_on_ejection<L>        ej;
	stop_on_close_encounter<L> ce;
	stop_on_crossing_orbit<L>  co;
};


/** Combination of stop_on_ejcetion and stop_on_close_encounter
 *
 */
template <class L> 
struct stop_on_ejection_or_close_encounter {
	struct params {
		typename stop_on_ejection<L>        ::params ej;
		typename stop_on_close_encounter<L> ::params ce;
		
		params(const config& cfg)
			:ej(cfg), ce(cfg) {}
	};
	
	GPUAPI stop_on_ejection_or_close_encounter
		(const params& p,ensemble::SystemRef& s,L& l)
		: ej(p.ej,s,l), ce(p.ce,s,l) 	{}
	
  //	GPUAPI void operator () () 
  //        {  ej(); ce(); }

	GPUAPI void operator () (const int thread_in_system) 
          {	ej(thread_in_system); ce(thread_in_system); 	}
	
private:
	stop_on_ejection<L>        ej;
	stop_on_close_encounter<L> ce;
};

  } } // end namespace monitors :: swarm



