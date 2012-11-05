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

/*! \file stop_on_any_large_distance.hpp
 *   \brief Defines and implements the monitor \ref swarm::monitors::stop_on_any_large_distance
 *          that logs the body that is separated from origin and other bodies by at least "rmax".        
 *
 *
 *  *EXPERIMENTAL*: This class is not thoroughly tested.
 * 
 */

#pragma once

#include <limits>

namespace swarm {
  namespace monitors {

/*! Parameters for stop_on_any_large_distance monitor
 * deactivate_on_close_encounter (bool): 
 * log_on_close_encounter (bool): 
 * verbose_on_close_encounter (bool): 
 * rmax (real): minimum distance between bodies to trigger
 *
 * \ingroup monitors_param
 *
 *  *EXPERIMENTAL*: This class is not thoroughly tested.
 * 
 */ 
struct stop_on_any_large_distance_params {
	double rmax;
        bool deactivate_on, log_on, verbose_on;
  /*! \param cfg Configuration Paramaters
   */
	stop_on_any_large_distance_params(const config &cfg)
	{
		rmax = cfg.optional("rmax",std::numeric_limits<float>::max());
		deactivate_on = cfg.optional("deactivate_on_any_large_distance",false);
		log_on = cfg.optional("log_on_any_large_distance",false);
		verbose_on = cfg.optional("verbose_on_any_large_distance",false);
	}


};

/** Simple monitor that logs when any one body is separated from 
 *  *EXPERIMENTAL*: This class is not thoroughly tested.
 *  \ingroup experimental
 *  both the origin and every other body by a distance of at least "rmax" 
 *  Optionally signals if "stop on rmax" is true
 *  This monitor may be useful in simple scattering experiments.
 *  WARNING: The code so that there is run time option of signaling to stop is untested.  
 *  TODO: Once it works, may want to port to other monitors.
 *  \ingroup monitors
 */
template<class log_t>
class stop_on_any_large_distance {
	public:
	typedef stop_on_any_large_distance_params params;

	private:
	params _params;
        bool need_full_test, condition_met;

	ensemble::SystemRef& _sys;
	log_t& _log;

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

        //! Check to see if need_full_test
	GPUAPI bool pass_one (int thread_in_system) 
          {
	    need_full_test = false; 
	    condition_met = false;
	    if(is_any_on()&&(thread_in_system==0))
	      {
		bool is_any_body_far_from_origin = false;
		for(int b = 0 ; b < _sys.nbod(); b ++ )
		  {
		    if(_sys.distance_to_origin_squared(b) > _params.rmax * _params.rmax )
		      is_any_body_far_from_origin = true;
		  }
		if(!is_any_body_far_from_origin) break;
		bool need_to_log = false;
		for(int b = 0 ; b < _sys.nbod(); b ++ )
		  {
		    if(_sys.distance_to_origin_squared(b) >= _params.rmax * _params.rmax ) 
		      {
			bool is_far_from_every_body = true;
			for(int bb = 0 ; bb < _sys.nbod(); bb ++ )
			  if(b != bb) 
			    {
			      double r2 = _sys.distance_squared_between(b,bb);
			      if(r2 < _params.rmax*_params.rmax )
				is_far_from_every_body = false;
			    }
			if(is_far_from_every_body) 
			  {
			    condition_met = true;
			    if( is_log_on() )
			      need_full_test = true;
			    if(is_verbose_on() )
			      lprintf(_log, "Distance from all bodies exceeds rmax: _sys=%d, bod=%d, T=%lg r=%lg rmax=%lg.\n", _sys.number(), b, _sys.time() , sqrt(r2), _params.rmax);
			  }
		      } // if far from origin
		  } // for bodies
	      }
	    return need_full_test;
	  }
	    
        //! Check the system state
	GPUAPI int pass_two (int thread_in_system) 
          {
	    if(is_condition_met() && is_deactivate_on() &&(thread_in_system==0) )
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
  GPUAPI void operator () (int thread_in_system) {
	  if(!is_any_on()) return;

		if(thread_in_system==0)
		  {

		bool is_any_body_far_from_origin = false;
		for(int b = 0 ; b < _sys.nbod(); b ++ ){
			if(_sys.distance_to_origin_squared(b) > _params.rmax * _params.rmax )
				is_any_body_far_from_origin = true;
		}
		if(!is_any_body_far_from_origin) return;
		bool need_to_log = false;
		for(int b = 0 ; b < _sys.nbod(); b ++ ){
			if(_sys.distance_to_origin_squared(b) >= _params.rmax * _params.rmax ) {
				bool is_far_from_every_body = true;
				for(int bb = 0 ; bb < _sys.nbod(); bb ++ )
					if(b != bb) {
						double r2 = _sys.distance_squared_between(b,bb);
						if(r2 < _params.rmax*_params.rmax )
							is_far_from_every_body = false;
					}
				if(is_far_from_every_body) {
				  if(is_verbose_on() )
					lprintf(_log, "Distance from all bodies exceeds rmax: _sys=%d, bod=%d, T=%lg r=%lg rmax=%lg.\n", _sys.number(), b, _sys.time() , sqrt(r2), _params.rmax);
				  if(is_log_on() )
				    need_to_log = true;

				  if(is_deactivate_on() )
					a_sys.set_disabled();
				}
			}
		}
		if(need_to_log)
		   log::system(_log, _sys);		  
		  }
  }
#endif

	GPUAPI stop_on_any_large_distance(const params& p,ensemble::SystemRef& s,log_t& l)
		:_params(p),_sys(s),_log(l)
#if 0
		,_counter(0)
#endif
                {}
	
};

}
}
