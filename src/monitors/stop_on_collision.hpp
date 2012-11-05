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

/*! \file stop_on_collision.hpp
 *   \brief Defines and implements the monitor \ref swarm::monitors::stop_on_collision
 *          that detects physical collisions. 
 *
 *  *EXPERIMENTAL*: This class does not pass any tests. 
 * 
 *
 */

#pragma once

namespace swarm {  namespace monitors {

/** Parameters for stop_on_collision monitor
 * deactivate_on_collision (bool): 
 * log_on_collision (bool): 
 * verbose_on_collision (bool): 
 * collision_distance_to_origin (real): default distance or collision if individual radii not avaliable
 *
 * \ingroup monitors_param
 */ 
struct stop_on_collision_param {
	double dmin_squared;
  bool deactivate_on, log_on, verbose_on;

  /*! \param cfg Configuration Paramaters
   */
	stop_on_collision_param(const config &cfg)
	{
	  dmin_squared = cfg.optional("collision_distance_to_origin",0.);
	  dmin_squared *= dmin_squared;

	  deactivate_on = cfg.optional("deactivate_on_collision",false);
	  log_on = cfg.optional("log_on_collision",false);
	  verbose_on = cfg.optional("verbose_on_collision",false);
	}
};

/** Simple monitor to detect physical collisions.  
 *  *EXPERIMENTAL*: This class is not thoroughly tested.
 *  \ingroup experimental
 *
 *  Signals and logs if current separation between any two bodies is less than "collision_distance_to_origin".
 *  WARNING: Does not interpolate between steps
 *
 *  \ingroup monitors
 */
template<class log_t>
class stop_on_collision {
	public:
	typedef stop_on_collision_param params;

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

        //! check for close encounters and return status of need_full_test
	GPUAPI bool pass_one (int thread_in_system) 
          {
	    need_full_test = false; 
	    condition_met = false;
	    if(is_any_on()&&(thread_in_system==0))
	      {
		// Chcek for close encounters
		for(int b = 1; b < _sys.nbod(); b++)
		  for(int d = 0; d < b; d++)
		    condition_met = condition_met || check_close_encounters(b,d); 
		if(condition_met && is_log_on() )
		    {  need_full_test = true;  }			

	      }
	    return need_full_test;
	  }
	    
        //! Check if deactivate the system and return the status
	GPUAPI int pass_two (int thread_in_system) 
          {
	    if (need_to_deactivate() && (thread_in_system()==0) )
	      {  _sys.set_disabled();  }
	    return _sys.state();
	  }


#if 0 // still under development
  GPUAPI bool check_close_encounter_possible(const int& i, const int& j, double dt){

    //	  double hill_distance_to_origin_sq_upper_limit = pow(((_sys.mass(i)+_sys.mass(j))/(3.0*_sys.mass(0))),2.0/3.0)*(std::max(_sys.distance_to_origin_squared(i),_sys.distance_to_origin_squared(j)));
	  double target_distance_to_origin_sq = (NUM_ATTRIBUTES>=1) ? attribute(i)*attribute(i) : _params.dmin_squared;
	  double vesc_sq = 2.0*_sys.mass(0)/std::min(_sys.distance_to_origin(i),_sys.distance_to_origin(j));
	  double d_squared = _sys.distance_squared_between(i,j);
	  bool close_encounter = (d_squared < vesc_sq*dt*dt + target_distance_to_origin_sq);
	  return close_encounter;
	}
#endif
  
#if 0 // still under development
  GPUAPI double min_encounter_distance(const int& i, const int& j){
    const Body& b1 = _body[i], & b2 = _body[j];
    double d = sqrt( square(b1[0].pos()-b2[0].pos())
		     + square(b1[1].pos()-b2[1].pos())
		     + square(b1[2].pos()-b2[2].pos()) );
    return d;
  }
#endif

        //! Check close encounters
	GPUAPI bool check_close_encounters(const int& i, const int& j){

		double d_squared = _sys.distance_squared_between(i,j);
		double target_distance_to_origin_sq = (NUM_ATTRIBUTES>=1) ? attribute(i)*attribute(i) : _params.dmin_squared;
		bool close_encounter = d_squared < target_distance_to_origin_sq;

		if( close_encounter )
		  if(is_verbose_on() )
			lprintf(_log, "Collision detected: "
					"sys=%d, T=%f j=%d i=%d  d=%lg.\n"
				, _sys.number(), _sys.time(), j, i,sqrt(d_squared));
		return close_encounter;
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
	GPUAPI void operator () (int thread_in_system) 
	  if(!is_any_on()) return;
	  bool stopit = false;

		if(thread_in_system==0)
		  {

		// Chcek for close encounters
		for(int b = 1; b < _sys.nbod(); b++)
			for(int d = 0; d < b; d++)
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

        //! default constructor for stop_on_collision
	GPUAPI stop_on_collision(const params& p,ensemble::SystemRef& s,log_t& l)
	    :_params(p),_sys(s),_log(l),_counter(0){}
	
};

} } // end namespace monitors :: swarm
