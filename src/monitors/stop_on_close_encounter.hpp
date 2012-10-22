/*************************************************************************
 * Copyright (C) 2011 by Eric Ford and the Swarm-NG Development Team     *
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

/*! \file stop_on_close_encounter.hpp
 *   \brief Defines and implements the monitor that signals and logs when the distance between
 *          any two bodies is less than "close_approach". 
 *
 */

#pragma once

#include <limits>

namespace swarm { namespace monitors {

/* Parameters for stop_on_close_encounter monitor
 * deactivate_on_close_encounter (bool): 
 * log_on_close_encounter (bool): 
 * verbose_on_close_encounter (bool): 
 * close_approach (real): maximum distance in Hill radii to trigger action
 *
 * \ingroup monitors_param
 */ 
struct stop_on_close_encounter_param {
	double dmin;
  bool deactivate_on, log_on, verbose_on;
  /*! \param cfg Configuration Paramaters
   */
	stop_on_close_encounter_param(const config &cfg)
	{
		dmin = cfg.optional("close_approach",0.0);
		deactivate_on = cfg.optional("deactivate_on_close_encounter",false);
		log_on = cfg.optional("log_on_close_encounter",false);
		verbose_on = cfg.optional("verbose_on_close_encounter",false);
	}
};

/** Simple monitor to detect close encounters.
 *  Signals and logs if current separation between any two bodies (measured in mutual Hill radii) is less than "close_approach".
 *  WARNING: Does not interpolate between steps
 *
 *  \ingroup monitors
 *  \ingroup monitors_for_planetary_systems
 */
template<class log_t>
class stop_on_close_encounter {
	public:
	typedef stop_on_close_encounter_param params;

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

	GPUAPI bool pass_one (int thread_in_system) 
          {
	    need_full_test = false; 
	    condition_met = false;
	    if(is_any_on()&&(thread_in_system==0))
	      {
		// Chcek for close encounters
		for(int b = 2; b < _sys.nbod(); b++)
		  for(int d = 1; d < b; d++)
		    condition_met = condition_met || check_close_encounters(b,d); 
		if( condition_met && is_log_on() )
		  {  need_full_test = true;  }

	      }
	    return need_full_test;
	  }
	    

	GPUAPI int pass_two (int thread_in_system) 
          {
	    if (need_to_deactivate() && (thread_in_system==0) )
	      {  _sys.set_disabled();  }
	    return _sys.state();
	  }


	GPUAPI bool check_close_encounters(const int& i, const int& j){

		double d = _sys.distance_between(i,j);
		double _GM = _sys[0].mass();  // remove _ if ok to keep
		//		double rH = pow((_sys[i].mass()+_sys[j].mass())/(3.*_GM),1./3.);
		//		bool close_encounter = d < _p.dmin * rH;
		double a = 0.5*(_sys[i].radius()+_sys[i].radius());
		double rH3 = (_sys[i].mass()+_sys[j].mass())/(3.*_GM)*a*a*a;
		bool close_encounter = d*d*d < _params.dmin*_params.dmin*_params.dmin * rH3;

		if( close_encounter )
		  if(is_verbose_on() )
			lprintf(_log, "Close apporach detected: "
					"sys=%d, T=%f j=%d i=%d  d=%lg.\n"
					, _sys.number(), _sys.time(), j, i,d);

		return close_encounter;
	}

        GPUAPI void operator () (const int thread_in_system) 
          { 
	    pass_one(thread_in_system);
	    pass_two(thread_in_system);
	    if(need_to_log_system() && (thread_in_system==0) )
	      log::system(_log, _sys);
	  }

#if 0
  //	GPUAPI void operator () ()  
	GPUAPI void operator () (int thread_in_system) 
  {
	  if(!is_any_on()) return;
		bool stopit = false;

		if(thread_in_system==0)
		  {
		// Chcek for close encounters
		for(int b = 2; b < _sys.nbod(); b++)
			for(int d = 1; d < b; d++)
				stopit = stopit || check_close_encounters(b,d); 

		if(stopit) {
		  if(is_log_on())
			log::system(_log, _sys);
		  if(is_deactivate_on())
			_sys.set_disabled();
		}
		  }
  }
#endif

	GPUAPI stop_on_close_encounter(const params& p,ensemble::SystemRef& s,log_t& l)
	    :_params(p),_sys(s),_log(l){}
	
};

} } // end namespace monitors :: swarm
