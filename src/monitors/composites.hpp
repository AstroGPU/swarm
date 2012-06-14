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
 * \ingroup monitors
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

#if 0
	GPUAPI void operator () (const int thread_in_system) 
         {  ej(thread_in_system); ce(thread_in_system); co(thread_in_system);  }
#endif

        GPUAPI void operator () (const int thread_in_system) 
          { 
	    pass_one(thread_in_system);
	    pass_two(thread_in_system);
	    if(need_to_log_system() && (thread_in_system==0) )
	      ej.log_system();
	  }
	
        GPUAPI bool is_deactivate_on() { return ej.is_deactivate_on() || ce.is_deactivate_on() || co.is_deactivate_on(); }
        GPUAPI bool is_log_on() { return ej.is_log_on() || ce.is_log_on() || co.is_log_on(); }
        GPUAPI bool is_verbose_on() { return ej.is_verbose_on() || ce.is_verbose_on() || co.is_verbose_on(); };
        GPUAPI bool is_any_on() { return is_deactivate_on() || is_log_on() || is_verbose_on() ; }

        GPUAPI bool is_condition_met () { return ej.is_condition_met() || ce.is_condition_met() || co.is_condition_met(); }
        GPUAPI bool need_to_log_system () 
          { return (ej.need_to_log_system() || ce.need_to_log_system() || co.need_to_log_system() ); }
        GPUAPI bool need_to_deactivate () 
           { return (ej.need_to_deactivate() || ce.need_to_deactivate()  || co.need_to_deactivate() ); }

	GPUAPI bool pass_one (int thread_in_system) 
          {  return ej.pass_one(thread_in_system) || ce.pass_one(thread_in_system) || co.pass_one(thread_in_system);  }
	    

	GPUAPI int pass_two (int thread_in_system) 
          {  
	    int s1 = ej.pass_two(thread_in_system);
	    int s2 = ce.pass_two(thread_in_system); 
	    int s3 = co.pass_two(thread_in_system); 
	    if((s1==0)&&(s1==0)&&(s3==0)) return 0;
	    else 
	      {
		int ret = s1;
		if((s1<0)||(s2<0)||(s3<0)) 
		  {
		    if(s2<ret) ret = s2;
		    if(s3<ret) ret = s3;
		  }
		else 
		  {
		    if(s2>ret) ret = s2;
		    if(s3>ret) ret = s3;
		  }
		return ret;
	      }
	  }

private:
	stop_on_ejection<L>        ej;
	stop_on_close_encounter<L> ce;
	stop_on_crossing_orbit<L>  co;
};


/** Combination of stop_on_ejcetion and stop_on_close_encounter
 *
 *
 * \ingroup monitors
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

#if 0
	GPUAPI void operator () (const int thread_in_system) 
          {	ej(thread_in_system); ce(thread_in_system); 	}
#endif
        GPUAPI void operator () (const int thread_in_system) 
          { 
	    pass_one(thread_in_system);
	    pass_two(thread_in_system);
	    if(need_to_log_system() && (thread_in_system==0) )
	      ej.log_system();
	  }
	
        GPUAPI bool is_deactivate_on() { return ej.is_deactivate_on() || ce.is_deactivate_on(); }
        GPUAPI bool is_log_on() { return ej.is_log_on() || ce.is_log_on(); }
        GPUAPI bool is_verbose_on() { return ej.is_verbose_on() || ce.is_verbose_on(); };
        GPUAPI bool is_any_on() { return is_deactivate_on() || is_log_on() || is_verbose_on() ; }

        GPUAPI bool is_condition_met () { return ej.is_condition_met() || ce.is_condition_met(); }
        GPUAPI bool need_to_log_system () 
          { return (ej.need_to_log_system() || ce.need_to_log_system() ); }
        GPUAPI bool need_to_deactivate () 
           { return (ej.need_to_deactivate() || ce.need_to_deactivate()); }


	GPUAPI bool pass_one (int thread_in_system) 
          {  return ej.pass_one(thread_in_system) || ce.pass_one(thread_in_system);  }
	    

	GPUAPI int pass_two (int thread_in_system) 
          {  
	    int s1 = ej.pass_two(thread_in_system);
	    int s2 = ce.pass_two(thread_in_system); 
	    if((s1==0)&&(s1==0)) return 0;
	    else 
	      {
		int ret = s1;
		if((s1<0)||(s2<0)) 
		  {  if(s2<ret) ret = s2;  }
		else 
		  {  if(s2>ret) ret = s2;  }
		return ret;
	      }
	  }

private:
	stop_on_ejection<L>        ej;
	stop_on_close_encounter<L> ce;
};

  } } // end namespace monitors :: swarm



