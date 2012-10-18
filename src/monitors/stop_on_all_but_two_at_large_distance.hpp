/*************************************************************************
 * Copyright (C) 2011 by Eric Ford and the Swarm-NG Development Team  *
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

/*! \file stop_on_all_but_two_at_large_distance.hpp
 *   \brief Defines and implements the monitor that signals and logs when no more 
 *          than two bodies are within a distance "rmax" of origin or another body.
 *
 */

#pragma once

#include <limits>

namespace swarm {
  namespace monitors {


/* Parameters for stop_on_all_but_two_at_large_distance monitor
 * deactivate_on_close_encounter (bool): 
 * log_on_close_encounter (bool): 
 * verbose_on_close_encounter (bool): 
 * rmax (real): minimum distance between bodies to be considered isolated
 *
 * \ingroup monitors_param
 */ 
struct stop_on_all_but_two_at_large_distance_params {
	double rmax;
        bool deactivate_on, log_on, verbose_on;
  /*! \param cfg Configuration Paramaters
   */
	stop_on_all_but_two_at_large_distance_params(const config &cfg)
	{
		double rmax = cfg.optional("rmax",std::numeric_limits<float>::max());
		deactivate_on = cfg.optional("deactivate_on_all_but_two_at_large_distance",false);
		log_on = cfg.optional("log_on_all_but_two_at_large_distance",false);
		verbose_on = cfg.optional("verbose_on_all_but_two_at_large_distance",false)
	}
};

/** Simple monitor that signals and logs when no more than two bodies are within a distance "rmax" of origin or another body.  Need to test this monitor.
 *  \ingroup monitors
 */
template<class log_t>
class stop_on_all_but_two_at_large_distance {
	public:
	typedef stop_on_all_but_two_at_large_distance_params params;

	private:
	params _params;
        bool need_full_test, condition_met;

	ensemble::SystemRef& _sys;
	log_t& _log;
	int _counter;

	public:
		template<class T>
		static GENERIC int thread_per_system(T compile_time_param){
			return 1;
		}

		template<class T>
		static GENERIC int shmem_per_system(T compile_time_param) {
			 return 0;
		}
        GPUAPI bool is_deactivate_on() { return _params.deactivate_on; };
        GPUAPI bool is_log_on() { return _params.log_on; };
        GPUAPI bool is_verbose_on() { return _params.verbose_on; };
        GPUAPI bool is_any_on() { return is_deactivate_on() || is_log_on() || is_verbose_on() ; }
        GPUAPI bool is_condition_met () { return ( condition_met ); }
        GPUAPI bool need_to_log_system () 
          { return (is_log_on() && is_condition_met() ); }
        GPUAPI bool need_to_deactivate () 
          { return ( is_deactivate_on() && is_condition_met() ); }

        GPUAPI void log_system()  {  log::system(_log, _sys);  }

	GPUAPI bool pass_one (int thread_in_system) 
          {
	    need_full_test = false; 
	    condition_met = false;
	    if(is_any_on()&&(thread_in_system==0))
	      {
		// Check for distance from origin
		int num_body_near_origin = 0, id1 = -1, id2 = -2;
		for(int b = 0 ; b < _sys.nbod(); b ++ ){
		  if(_sys.radius_squared(b) <= _params.rmax*_params.rmax ) // WARNING: Confusing function name
		    {
		      if(num_body_near_origin==0) id1 = b;
		      if(num_body_near_origin==1) id2 = b;
		      num_body_near_origin++;
		    }
		}
		if( num_body_near_origin <= 2 )
		  {  need_full_test = true;  }
	      }
	    return need_full_test;
	  }
	    

	GPUAPI int pass_two (int thread_in_system) 
          {
	    int new_state = _sys.state();
	    if(is_any_on() && need_full_test) // implies thread_in_system==0
	      {
		// Check for distance from other bodies
		int num_body_far_from_all = 0;
		for(int b = 0 ; b < _sys.nbod(); b ++ )
		  {
		    if(_sys.radius_squared(b) <= _params.rmax*_params.rmax ) continue; // WARNING: Confusing function name
		    bool is_far_from_every_body = true;
		    for(int bb = 0 ; bb < _sys.nbod(); bb ++ ){		
		      if(b == bb) continue;
		      if(_sys.distance_squared_between(b,bb) < _params.rmax*_params.rmax )
			{ is_far_from_every_body = false;  break; }
		    }
		    if(is_far_from_every_body) 
		      {   num_body_far_from_all++;	  }
		  }
		if(num_body_far_from_all+2>=_sys.nbod())
		  {
		    condition_met = true;
		    if(is_verbose_on())
		      lprintf(_log, "No more than two bodies are within rmax of other bodies: _sys=%d, bod1=%d, bod2=%d, T=%lg r=%lg rmax=%lg.\n"
			      , _sys.number(), id1, id2, _sys.time() , r, _params.rmax);
		    if(is_deactivate_on() && (_sys.state()>=0) )
		      {  _sys.set_disabled();  new_state = _sys.state(); }
		  }
	      }
	    return new_state;
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
  GPUAPI void operator () (int thread_in_system) {
	  if(!is_any_on()) return;

		if(thread_in_system==0)
		  {
		// Check for distance from origin
		int num_body_near_origin = 0, id1 = -1, id2 = -2;
		for(int b = 0 ; b < _sys.nbod(); b ++ ){
			if(_sys.radius_squared(b) <= _params.rmax*_params.rmax ) // WARNING: Confusing function name
				{
				if(num_body_near_origin==0) id1 = b;
				if(num_body_near_origin==1) id2 = b;
				 num_body_near_origin++;
				}
		}
		if( num_body_near_origin > 2 )
		  return;
 
		// Check for distance from other bodies
		int num_body_far_from_all = 0;
		for(int b = 0 ; b < _sys.nbod(); b ++ ){
			if(_sys.radius_squared(b) <= _params.rmax*_params.rmax ) continue; // WARNING: Confusing function name
			bool is_far_from_every_body = true;
			for(int bb = 0 ; bb < _sys.nbod(); bb ++ ){		
			   if(b == bb) continue;
			   if(_sys.distance_squared_between(b,bb) < _params.rmax*_params.rmax )
					{ is_far_from_every_body = false;  break; }
				}
			if(is_far_from_every_body) 
				{
				num_body_far_from_all++;
				}
		}
		if(num_body_far_from_all+2>=_sys.nbod())
			{
			  if(is_verbose_on())
			lprintf(_log, "No more than two bodies are within rmax of other bodies: _sys=%d, bod1=%d, bod2=%d, T=%lg r=%lg rmax=%lg.\n"
				, _sys.number(), id1, id2, _sys.time() , r, _params.rmax);
			  if(is_log_on())
			log::system(_log, _sys);
			  if(is_deactivate_on())
			  _sys.set_disabled();
			}
		  }
  }
#endif

	GPUAPI stop_on_all_but_two_at_large_distance(const params& p,ensemble::SystemRef& s,log_t& l)
		:_params(p),_sys(s),_log(l),_counter(0){}
	
};

  }
}

