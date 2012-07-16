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

#define USE_WORKS 1

namespace swarm {
  namespace monitors {

/* Parameters for log_transit monitor
 * log_transit_tol (real): desired precision of transit times (zero makes inactive)
 * 
 * \ingroup monitors_param
 */ 
struct log_transit_params {
        double tol, dt;
        int nbod, num_max_iter;
        bool log_on_transit, log_on_occultation;

	log_transit_params(const config &cfg)
	{
	  nbod = cfg.require<int>("nbod"); 
	  num_max_iter = cfg.optional("num_max_transit_iter", 2); 
	  tol = cfg.optional("log_transit_tol", 2.e-8); // 0.1 seconds for G=Msol=AU=1
	  // WARNING assumes Hermite Fixed
	  dt = cfg.require("time_step", 0.0);  // use hermite's timestep
	  log_on_transit = cfg.optional("log_transits", true); 
	  log_on_occultation = cfg.optional("log_occultations", false); 
	  
	}
};

/** Monitor that logs (the estimated transit time, the estimated minimum impact parameter in units of stellar radii, and the estimated velocity projected onto the planet of the sky in units of stellar radii per time) at times near a transit
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
	  }
	    

  GPUAPI double square(const double x) { return x*x; }

  template<int nbod>
  GPUAPI void calc_transit_time(const int& thread_in_system, const int& i,const int& j, const double& dt, double dx[2], double dv[2], const double& b2begin, double& db2dt, const double& pos_step_end, const double& vel_step_end, double& dt_min_b2,double& b,double & vproj)
  {
    extern __shared__ char shared_mem[];
    typedef typename swarm::gpu::bppt::Gravitation<nbod>::shared_data grav_t;
    typedef swarm::gpu::bppt::Gravitation<nbod> calcForces_t;
    calcForces_t calcForces(_sys,*( (grav_t*) (&shared_mem[(swarm::gpu::bppt::sysid_in_block()%SHMEM_CHUNK_SIZE)*sizeof(double)+(swarm::gpu::bppt::sysid_in_block()/SHMEM_CHUNK_SIZE)*SHMEM_CHUNK_SIZE*(nbod*(nbod-1)*3*sizeof(double))]) ) );

    const int bid = swarm::gpu::bppt::thread_body_idx(nbod);
    const int cid = swarm::gpu::bppt::thread_component_idx(nbod);
    double pos = pos_step_end, vel = vel_step_end;
    double acc, jerk;
#if USING_INTEGRATOR_THAT_OVER_WRITES_SHARED_DATA || 1
        calcForces(thread_in_system,bid,cid,pos,vel,acc,jerk);
	/*    __syncthreads();
    if(thread_in_system < calcForces_t::pair_count)
     calcForces.calc_pair(thread_in_system);
     __syncthreads(); */
#endif
    double acc_ix, acc_jx, acc_iy, acc_jy, jerk_ix, jerk_jx, jerk_iy, jerk_jy;
    calcForces.sum(i,0,acc_ix,jerk_ix);
    calcForces.sum(i,1,acc_iy,jerk_iy);
    calcForces.sum(j,0,acc_jx,jerk_jx);
    calcForces.sum(j,1,acc_jy,jerk_jy);
    double da[2] = { acc_jx - acc_ix,  acc_jy- acc_iy };
    double dj[2] = { jerk_jx-jerk_ix, jerk_jy-jerk_iy };
    double d2b2dt2 = 2.*(dv[0]*dv[0]+dv[1]*dv[1]+dx[0]*da[0]+dx[1]*da[1]);
    double d3b2dt3 = 6.*(dv[0]*da[0]+dv[1]*da[1])+2.*(dx[0]*dj[0]+dx[1]*dj[1]);
    double dt_try = -db2dt/d2b2dt2;  // find vertex of parabola
    dt_try = -db2dt/(d2b2dt2+0.5*dt_try*d3b2dt3);
    // if((dt_try<0.) || (dt_try>dt) )             // ignore if already past or vertex past next step
    if((dt_try<=-0.5*dt) || (dt_try>0.5*dt) ) // ignore if already past or vertex past next step
      { dt_min_b2 = dt_try; b=-1.; vproj = -1.; return;  }
    double dt_cum = dt_try;
    double b2;
    for(int iter=0;iter<_params.num_max_iter;++iter)
      {
	// take trial step towards transit mid-point
	//	pos = pos_step_end + dt_try*(vel+(dt_try*0.5)*(acc+(dt_try/3.0)*jerk));
	//	vel = vel_step_end + dt_try*(acc+(dt_try*0.5)*jerk);
	pos = pos + dt_try*(vel+(dt_try*0.5)*(acc+(dt_try/3.0)*jerk));
	vel = vel + dt_try*(acc+(dt_try*0.5)*jerk);
	calcForces(thread_in_system,bid,cid,pos,vel,acc,jerk);
	dx[0] = _sys[j][0].pos()-_sys[i][0].pos();
	dx[1] = _sys[j][1].pos()-_sys[i][1].pos();
	dv[0] = _sys[j][0].vel()-_sys[i][0].vel();
	dv[1] = _sys[j][1].vel()-_sys[i][1].vel();
	b2 = dx[0]*dx[0]+dx[1]*dx[1];
	db2dt = 2.*(dx[0]*dv[0]+dx[1]*dv[1]);
	calcForces.sum(i,0,acc_ix,jerk_ix);
	calcForces.sum(i,1,acc_iy,jerk_iy);
	calcForces.sum(j,0,acc_jx,jerk_jx);
	calcForces.sum(j,1,acc_jy,jerk_jy);
	da[0] = acc_jx - acc_ix; da[1] = acc_jy- acc_iy;
	dj[0] = jerk_jx-jerk_ix, dj[1] = jerk_jy-jerk_iy;
	d2b2dt2 = 2.*(dv[0]*dv[0]+dv[1]*dv[1]+dx[0]*da[0]+dx[1]*da[1]);
	d3b2dt3 = 6.*(dv[0]*da[0]+dv[1]*da[1])+2.*(dx[0]*dj[0]+dx[1]*dj[1]);
	double dt_try2 = -db2dt/d2b2dt2;  // find vertex of parabola
	dt_try2 =  -db2dt/(d2b2dt2+0.5*dt_try2*d3b2dt3);
	//	dt_try += dt_try2;
	dt_cum += dt_try2;
	dt_try = dt_try2;
	if(dt_try2>_params.tol) continue;
	b2 += dt_try2*(db2dt+0.5*dt_try2*(d2b2dt2+dt_try2*d3b2dt3/3.));
	dv[0] += dt_try2*(da[0]+0.5*dt_try2*dj[0]);
	dv[1] += dt_try2*(da[1]+0.5*dt_try2*dj[1]);
	if(dt_try2<=_params.tol) break;
      }
    dt_try = dt_cum;
    // if((dt_try<0.) || (dt_try>dt) )             // ignore if already past or vertex past next step
    if((dt_try<=-0.5*dt) || (dt_try>0.5*dt) ) // ignore if already past or vertex past next step
      { dt_min_b2 = dt_try; b=-1.; vproj = -1.; return;  }
    
    dt_min_b2 = 0.;
    double min_b2 = b2begin;
    if(b2<min_b2) { min_b2 = b2; dt_min_b2 = dt_try; }
    b = (min_b2>=0) ? sqrt(min_b2) : -sqrt(-min_b2);
    vproj = sqrt(dv[0]*dv[0]+dv[1]*dv[1]);
  }


  // save one based on simple extrapolation, while working on one that itterates
  template<int nbod>
  GPUAPI void calc_transit_time_works(const int& i,const int& j, const double& dt, const double dx[2], double dv[2], const double& b2begin, double& db2dt,double& dt_min_b2,double& b,double & vproj)
  {
    extern __shared__ char shared_mem[];
    typedef typename swarm::gpu::bppt::Gravitation<nbod>::shared_data grav_t;
    swarm::gpu::bppt::Gravitation<nbod> calcForces(_sys,*( (grav_t*) (&shared_mem[(swarm::gpu::bppt::sysid_in_block()%SHMEM_CHUNK_SIZE)*sizeof(double)+(swarm::gpu::bppt::sysid_in_block()/SHMEM_CHUNK_SIZE)*SHMEM_CHUNK_SIZE*(nbod*(nbod-1)*3*sizeof(double))]) ) );

    // if integrator could have overwritten data in shared memory, need new call to calcforces
    double acc_ix, acc_jx, acc_iy, acc_jy, jerk_ix, jerk_jx, jerk_iy, jerk_jy;
    calcForces.sum(i,0,acc_ix,jerk_ix);
    calcForces.sum(i,1,acc_iy,jerk_iy);
    calcForces.sum(j,0,acc_jx,jerk_jx);
    calcForces.sum(j,1,acc_jy,jerk_jy);
    double da[2] = { acc_jx - acc_ix,  acc_jy- acc_iy };
    double dj[2] = { jerk_jx-jerk_ix, jerk_jy-jerk_iy };
    double d2b2dt2 = 2.*(dv[0]*dv[0]+dv[1]*dv[1]+dx[0]*da[0]+dx[1]*da[1]);
    double d3b2dt3 = 6.*(dv[0]*da[0]+dv[1]*da[1])+2.*(dx[0]*dj[0]+dx[1]*dj[1]);
    double dt_try = -db2dt/d2b2dt2;  // find vertex of parabola
    dt_try = -db2dt/(d2b2dt2+0.5*dt_try*d3b2dt3);
    //    if((dt_try<0.) || (dt_try>dt) )                // ignore if already past or vertex past next step
    if((dt_try<-0.5*dt) || (dt_try>0.5*dt) ) // ignore if already past or vertex past next step
      { dt_min_b2 = dt_try; b=-1.; vproj = -1.; return;  }

    dt_min_b2 = 0.;
    double min_b2 = b2begin;
    double b2 = b2begin + dt_try*(db2dt + dt_try*0.5*(d2b2dt2+dt_try*d3b2dt3/3.));
    //dx[0] += dt_min_b2*(dv[0]+dt_min_b2*0.5*(da[0]+dt_min_b2*dj[0]/3.));
    //dx[1] += dt_min_b2*(dv[1]+dt_min_b2*0.5*(da[1]+dt_min_b2*dj[1]/3.));
    //double b2 = sqrt(dx[0]*dx[0]+dx[1]*dx[1]);
    if(b2<min_b2) { min_b2 = b2; dt_min_b2 = dt_try; }
    //    if(min_b2>square((radi+radj)*(1.+_params.tol)))       return false;
    b = (min_b2>0) ? sqrt(min_b2) : -sqrt(-min_b2);
    dv[0] += dt_min_b2*(da[0]+0.5*dt_min_b2*dj[0]);
    dv[1] += dt_min_b2*(da[1]+0.5*dt_min_b2*dj[1]);
    vproj = sqrt(dv[0]*dv[0]+dv[1]*dv[1]);
  }


  GPUAPI bool check_in_transit(const int& thread_in_system, const int& i, const int& j, const double dt)
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

	    // store values from actual integration to return to after finding transit time
	    const int nbod = _params.nbod;
	    const int bid = swarm::gpu::bppt::thread_body_idx(nbod);
	    const int cid = swarm::gpu::bppt::thread_component_idx(nbod);
	    double pos_step_end, vel_step_end;

	    if( (bid < nbod) && (cid < 3) ) {
	      pos_step_end = _sys[bid][cid].pos();
	      vel_step_end = _sys[bid][cid].vel();
	    }

	    // WARNING: hard coded to assume Hermite Fixed, nbod=3 and Gravitation
	    double dt_min_b2, b, vproj;
	    if((nbod==3) && (nbod<MAX_NBODIES))
	       {
		 const int nbod_hardwired = 3;
#if USE_WORKS
		 calc_transit_time_works<nbod_hardwired> (i,j,dt,dx,dv,b2begin,db2dt,dt_min_b2,b,vproj);
#else
		 calc_transit_time<nbod_hardwired> (thread_in_system,i,j,dt,dx,dv,b2begin,db2dt,pos_step_end,vel_step_end,dt_min_b2,b,vproj);
#endif
	       }
#if 0
	     else if((nbod==4) && (nbod<MAX_NBODIES))
	       {
		 const int nbod_hardwired = 4;
#if USE_WORKS
		 calc_transit_time_works<nbod_hardwired> (i,j,dt,dx,dv,b2begin,db2dt,dt_min_b2,b,vproj);
#else
		 calc_transit_time<nbod_hardwired> (thread_in_system,i,j,dt,dx,dv,b2begin,db2dt,pos_step_end,vel_step_end,dt_min_b2,b,vproj);
#endif
	       }
	     else if(nbod==5)
	       {
		 const int nbod_hardwired = 5;
#if USE_WORKS
		 calc_transit_time_works<nbod_hardwired> (i,j,dt,dx,dv,b2begin,db2dt,dt_min_b2,b,vproj);
#else
		 calc_transit_time<nbod_hardwired> (thread_in_system,i,j,dt,dx,dv,b2begin,db2dt,pos_step_end,vel_step_end,dt_min_b2,b,vproj);
#endif
	       }
	     else if((nbod==6) && (nbod<MAX_NBODIES))
	       {
		 const int nbod_hardwired = 6;
#if USE_WORKS
		 calc_transit_time_works<nbod_hardwired> (i,j,dt,dx,dv,b2begin,db2dt,dt_min_b2,b,vproj);
#else
		 calc_transit_time<nbod_hardwired> (thread_in_system,i,j,dt,dx,dv,b2begin,db2dt,pos_step_end,vel_step_end,dt_min_b2,b,vproj);
#endif
	       }
	     else if((nbod==7) && (nbod<MAX_NBODIES))
	       {
		 const int nbod_hardwired = 7;
#if USE_WORKS
		 calc_transit_time_works<nbod_hardwired> (i,j,dt,dx,dv,b2begin,db2dt,dt_min_b2,b,vproj);
#else
		 calc_transit_time<nbod_hardwired> (thread_in_system,i,j,dt,dx,dv,b2begin,db2dt,pos_step_end,vel_step_end,dt_min_b2,b,vproj);
#endif
	       }
	     else if((nbod==8) && (nbod<MAX_NBODIES))
	       {
		 const int nbod_hardwired = 8;
#if USE_WORKS
		 calc_transit_time_works<nbod_hardwired> (i,j,dt,dx,dv,b2begin,db2dt,dt_min_b2,b,vproj);
#else
		 calc_transit_time<nbod_hardwired> (thread_in_system,i,j,dt,dx,dv,b2begin,db2dt,pos_step_end,vel_step_end,dt_min_b2,b,vproj);
#endif
	       }
#endif
	     if(vproj<0.) 
	       {  // ignore if already past or vertex past next step
		 if( (bid < nbod) && (cid < 3) ) {
		   _sys[bid][cid].pos() = pos_step_end;
		   _sys[bid][cid].vel() = vel_step_end;
		 }
		 return false;
	       }
	     if(thread_in_system==0)
	       {
		 if(_sys.num_body_attributes()>=1) 
		   {
		     if (depth>0.)  { b /= radi; vproj /= radi; }
		     else      { b /= radj; vproj /= radj; }
		   }
		 
		 int event_id = (depth>0.) ? log::EVT_TRANSIT : log::EVT_OCCULTATION;
		 log::event(_log,event_id,_sys.time()+dt_min_b2,_sys.id(),j,b,vproj);
	       }
	     // restore values from main integration	
	    if( (bid < nbod) && (cid < 3) ) {
	      _sys[bid][cid].pos() = pos_step_end;
	      _sys[bid][cid].vel() = vel_step_end;
	    }

	     return true;
	  }
  

	GPUAPI int pass_two (const int thread_in_system) 
          {   
            if(is_any_on())
              {
                // Chcek for close encounters
                for(int b = 1; b < _sys.nbod(); b++)
		  {
		    // use with _works version
#if USE_WORKS
		    		    if(thread_in_system==0) 
		    		    condition_met = condition_met || 
#endif
		    check_in_transit(thread_in_system,0,b,_params.dt);
      		  }
	      }
	    return condition_met;
	  }

   

	GPUAPI log_transit(const params& p,ensemble::SystemRef& s,log_t& l)
	  :_params(p),_sys(s),_log(l) {}
	
};

}


}
