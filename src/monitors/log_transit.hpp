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
#pragma once

#include "../swarm/gpu/gravitation.hpp"

namespace swarm {
  namespace monitors {

/* Parameters for log_transit monitor
 * log_transit_tol (real): desired precision of transit times (zero makes inactive)
 * 
 * \ingroup monitors_param
 */ 
struct log_transit_params {
        double tol, dt;
        bool log_on_transit, log_on_occultation;

	log_transit_params(const config &cfg)
	{
	  tol = cfg.optional("log_transit_tol", 2.e-8); // 0.1 seconds for G=Msol=AU=1
	  // WARNING assumes Hermite Fixed
	  dt = cfg.require("time_step", 0.0);  // use hermite's timestep
	  log_on_transit = cfg.optional("log_on_transit", true); 
	  log_on_occultation = cfg.optional("log_on_occultation", false); 
	  
	}
};

/** Monitor that logs what? (the estimated transit time, the estimated minimum impact parameter in units of stellar radii, and the estimated velocity projected onto the planet of the sky in units of stellar radii per time) at times near a transit
 *   *  Assumes integration results in increasing time.
 *   *  Currently, highly experimental
 *   *  Currently, hardwired for Hermite, Fixed, std. Gravitation, nbod=3
 *   *  
 * 
 *  \ingroup monitors
 *
 *   parameters have yet to be finalized
 */
template<class log_t>
class log_transit {
	public:
	typedef log_transit_params params;

	private:
	params _params;
        bool condition_met;

	ensemble::SystemRef& _sys;
  //	double _next_log_time;
	log_t& _log;

	public:
        GPUAPI bool is_deactivate_on() { return false; };
  //        GPUAPI bool is_log_on() { return _params.log_on; };
        GPUAPI bool is_log_on() { return _params.tol!=0.; };
        GPUAPI bool is_verbose_on() { return false; };
        GPUAPI bool is_any_on() { return is_deactivate_on() || is_log_on() || is_verbose_on() ; }
        GPUAPI bool is_condition_met () { return ( condition_met ); }
        GPUAPI bool need_to_log_system () 
          { return (is_log_on() && is_condition_met() ); }
        GPUAPI bool need_to_deactivate () 
          { return ( is_deactivate_on() && is_condition_met() ); }

        GPUAPI void log_system()  {  log::system(_log, _sys);  }
	
        GPUAPI void operator () (const int thread_in_system) 
          { 
	    pass_one(thread_in_system);
	    pass_two(thread_in_system);
	    if(need_to_log_system() && (thread_in_system==0) )
	      log_system();
	  }

	GPUAPI bool pass_one (const int thread_in_system) 
          {
            condition_met = false;
	    return true;
#if 0	    // Could put code to only do full test if close to transit
	    need_full_test = false; 
            if(is_any_on()&&(thread_in_system==0))
              {
                // Chcek for close encounters
                for(int b = 1; b < _sys.nbod(); b++)
		  condition_met = condition_met || check_near_transit(0,b); 
                if(is_condition_met() && is_log_on() )
                    {  need_full_test = true;  }                        

              }
            return need_full_test;
#endif
	  }
	    

  GPUAPI double square(const double x) { return x*x; }
  GPUAPI bool check_in_transit(const int& i, const int& j, const double dt)
          {
	    double depth = _sys[j][2].pos()-_sys[i][2].pos();
	    if(!_params.log_on_occultation && (depth < 0.)) return false;
	    if(!_params.log_on_transit && (depth > 0.)) return false; 

	    double dx[2] = { _sys[j][0].pos()-_sys[i][0].pos(), _sys[j][1].pos()-_sys[i][1].pos() };
	    double dv[2] = { _sys[j][0].vel()-_sys[i][0].vel(), _sys[j][1].vel()-_sys[i][1].vel() };
	    double b2begin = dx[0]*dx[0]+dx[1]*dx[1];
	    double db2dt = 2.*(dx[0]*dv[0]+dx[1]*dv[1]);

	    double radi = (_sys.num_body_attributes()>=1) ? _sys[i].attribute(0) : 0.;
	    double radj = (_sys.num_body_attributes()>=1) ? _sys[j].attribute(0) : 0.;
	    if( min(b2begin,b2begin+dt*db2dt)>square((radi+radj)*(1.+_params.tol)) )
	    if( min(b2begin-0.5*dt*db2dt,b2begin+0.5*dt*db2dt)>square((radi+radj)*(1.+_params.tol)) )
	      return false;
	      
	    // WARNING: hard coded to assume Hermite Fixed, nbod=3 and Gravitation
	     extern __shared__ char shared_mem[];
	     const int nbod = 3;
             typedef typename swarm::gpu::bppt::Gravitation<nbod>::shared_data grav_t;
	     // assert(nbod==3);
             swarm::gpu::bppt::Gravitation<nbod> calcForces(_sys,*( (grav_t*) (&shared_mem[(swarm::gpu::bppt::sysid_in_block()%SHMEM_CHUNK_SIZE)*sizeof(double)+(swarm::gpu::bppt::sysid_in_block()/SHMEM_CHUNK_SIZE)*SHMEM_CHUNK_SIZE*(nbod*(nbod-1)*3*sizeof(double))]) ) );

	    double da[2] = { calcForces.sum_acc(j,0)-calcForces.sum_acc(i,0), calcForces.sum_acc(j,1)-calcForces.sum_acc(i,1) };
	    double d2b2dt2 = 2.*(dv[0]*dv[0]+dv[1]*dv[1]+dx[0]*da[0]+dx[1]*da[1]);
	    double dt_min_b2 = 0., min_b2 = b2begin;
	    double dt_try = -db2dt/d2b2dt2;  // find vertex of parabola
	    //	    if((dt_try<0.) || (dt_try>dt) ) return false;  // ignore if already past or vertex past next step
	    if((dt_try<-0.5*dt) || (dt_try>0.5*dt) ) return false;  // ignore if already past or vertex past next step
#if 1  // simple extrapolation
	    //	    db2dt = 2.*(dx[0]*dv[0]+dx[1]*dv[1])+dt_try*(2.0*(dv[0]*dv[0]+dv[1]*dv[1])+dx[0]*da[0]+dx[1]*da[1])+dt_try*dt_try*(3.*(dv[0]*da[0]+dv[1]*da[1])+dx[0]*dj[0]+dx[1]*dj[1]);
	    //      d2b2dt2 = 
	    double b2 = b2begin + dt_try*(db2dt + dt_try*0.5*d2b2dt2);
	    //	    double db2dt += dt*((2.*(dv[0]*dv[0]+dv[1]*dv[1])+dx[0]*da[0])+dt*((2.*(dv[0]*da[0]+dv[1]*da[1]))+dx[0]*da[0]+dx[1]*da[1]));
	    if(b2<min_b2) { min_b2 = b2; dt_min_b2 = dt_try; }
	    if(min_b2>square((radi+radj)*(1.+_params.tol))) return false;
	    //	    dx[0] += dt_min_b2*(dv[0]+dt_min_b2*da[0]);
	    //	    dx[1] += dt_min_b2*(dv[1]+dt_min_b2*da[1]);
	    double b = (min_b2>0) ? sqrt(min_b2) : 0.;
	    dv[0] += dt_min_b2*da[0];
	    dv[1] += dt_min_b2*da[1];
	    double vproj = sqrt(dv[0]*dv[0]+dv[1]*dv[1]);
#else
	    // include a propagator in monitor?!? to itterate to get transit time
#endif
	    if(_sys.num_body_attributes()>=1) 
	      {
		if (depth>0.)  { b /= radi; vproj /= radi; }
		else      { b /= radj; vproj /= radj; }
	      }
	    
	    int event_id = (depth>0.) ? log::EVT_TRANSIT : log::EVT_OCCULTATION;
	    log::event(_log,event_id,_sys.time()+dt_min_b2,_sys.id(),j,b,vproj);
	    log_system();	    
	    return true;
	  }


	GPUAPI int pass_two (const int thread_in_system) 
          {   
            if(is_any_on())
              {
                // Chcek for close encounters
                for(int b = 1; b < _sys.nbod(); b++)
		  {
		    if(thread_in_system==0)
		      condition_met = condition_met || check_in_transit(0,b,_params.dt);
		  }
	      }
	    return condition_met;
	  }

   

	GPUAPI log_transit(const params& p,ensemble::SystemRef& s,log_t& l)
	  :_params(p),_sys(s),_log(l) {}
	
};

}


}
