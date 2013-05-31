/** Implementation of the algorithm for collision detection in Mecurial 
 *	by John E. Chambers which could be found at 
 *  https://github.com/smirik/mercury
 *************************************************************************
 * Copyright (C) 2011 by Thien Nguyen and the Swarm-NG Development Team  *
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
 */
 
#pragma once

#include <limits>
#include <cstdlib>

#define PI 3.14159265358979
#define THIRD 0.3333333333333333
namespace swarm {

/*!
 * @brief Namespace for monitors (i.e., stoppers & loggers) and their associated classes
 */
namespace monitors { 
 
  //! Structure for monitor_template_params
 struct mce_stat_params {
  int nbod;
  double rceh, rho, dmin, time_step;
  bool deactivate_on, log_on, verbose_on;
  /*! \param cfg Configuration Paramaters
   */
	mce_stat_params(const config &cfg)
	{
		nbod = cfg.require<int>("nbod");
		time_step = cfg.require("time_step", 0.0);
		dmin = cfg.optional("close_approach",0.0);
		rceh = cfg.optional("rceh",0.001);
		rho = cfg.optional("rho",5.52); // density of the earth or planet in general
		
		deactivate_on = cfg.optional("deactivate_on_close_encounter",false);
		log_on = cfg.optional("log_on_close_encounter",false);
		verbose_on = cfg.optional("verbose_on_close_encounter",false);
		
	} 
   
 };

/** Empty monitor to use as a template.  
 * Signal is always false.  Does not do any logging.
 * \ingroup monitors
 *
 */
template<class log_t>
class mce_stat {
	public:
	typedef mce_stat_params params;

	private:
	params _params;
	ensemble::SystemRef& _sys;
	log_t& _log;
	
	double* rsc; // shared memory
	
	int nbod;
	int pair_count;
	
	double timestep;
	
	/// Bounding box of each particle between two position x0 at time t0 and x1 at time t1
	double *BBox;
	
	/// Position and velocity of the previous time step
	double *X0, *V0;
	
	// hill*rceh; hill radii could be calculated from mass, position, velocity; rceh is initial input
	//double *rce;
	
	// Output
	//int stat[nbod]; // status (=0 => alive, <>0 => to be removed)
	
	public:

	  //! default monitor_template construct
	GPUAPI mce_stat(const params& p,ensemble::SystemRef& s,log_t& l, double* shared_mem)
		:_params(p),_sys(s),_log(l), rsc(shared_mem)
		{
		    nbod = _params.nbod;
		    timestep = _params.time_step;
		    pair_count = (nbod*(nbod-1))/2;
		    BBox =  rsc;
		    X0 = rsc + nbod*3;
		    V0 = X0 + nbod*3;
		    
		}
	
	template<class T>
	static GENERIC int thread_per_system(T compile_time_param){
	  return (T::n*(T::n-1))/2;
	}
		
	template<class T>
	static GENERIC int shmem_per_system(T compile_time_param) {
		  return 12*T::n*sizeof(double);
	}
	
        //! Provide these functions, so two monitors can be combined
        GPUAPI bool is_deactivate_on() { return false; }
        GPUAPI bool is_log_on() { return false; }
        GPUAPI bool is_verbose_on() { return false; }
        GPUAPI bool is_any_on() { return is_deactivate_on() || is_log_on() || is_verbose_on() ; }
        GPUAPI bool is_condition_met () { return false; }
        GPUAPI bool need_to_log_system () 
          { return (is_log_on() && is_condition_met() ); }
        GPUAPI bool need_to_deactivate () 
          { return ( is_deactivate_on() && is_condition_met() ); }

    GPUAPI void init(const int thread_in_system)
    {
      if (thread_in_system < nbod)
      {
	_sys[thread_in_system].get(X0[3*thread_in_system],
				   X0[3*thread_in_system+1],
    				   X0[3*thread_in_system+2],
  				   V0[3*thread_in_system], 
    				   V0[3*thread_in_system+1],
				   V0[3*thread_in_system+2]);
	if (thread_in_system > 0)
	{
	  double gm = _sys[0].mass() + _sys[thread_in_system].mass();
	  double r = sqrt(X0[3*thread_in_system]*X0[3*thread_in_system] + X0[3*thread_in_system+1]*X0[3*thread_in_system+1] + X0[3*thread_in_system+2]*X0[3*thread_in_system+2]);
	  double v2 = V0[3*thread_in_system]*V0[3*thread_in_system] + V0[3*thread_in_system+1]*V0[3*thread_in_system+1] + V0[3*thread_in_system+2]*V0[3*thread_in_system+2];
	  double a = gm*r/(2*gm - r*v2);
	  if ( a <=0 ) a = r;
	  double hill = a*pow(THIRD*_sys[thread_in_system].mass()/_sys[0].mass(), THIRD);
	  _sys[thread_in_system].attribute(0) = hill*_params.rceh; // rce
	  _sys[thread_in_system].attribute(1) = hill/a * pow(2.25 * _sys[0].mass()/(PI*_params.rho),THIRD); // rphys
	}
	else
	{
	  lprintf(_log, "Initializing mce_stat monitor !\n");
	}
      }
    }

    GPUAPI void operator () (const int thread_in_system) 
    { 
    	
	    pass_one(thread_in_system);
	    
	    __syncthreads();
	    
	    pass_two(thread_in_system);
	    
	    __syncthreads();
	    
	    if(need_to_log_system() && (thread_in_system==0) )
	      log::system(_log, _sys);
	      
	    // store the current position and velocity for the next step
	    _sys[thread_in_system].get(X0[3*thread_in_system],
    				   X0[3*thread_in_system+1],
    				   X0[3*thread_in_system+2],
  				   V0[3*thread_in_system],
    				   V0[3*thread_in_system+1],
				   V0[3*thread_in_system+2]);
	
    	
    }

	private:
		
        /// This is mce_box() routine in Mercury,
        /// Calculate bounding box of the trajectory of each planet between time steps t0 and t1
	GPUAPI void pass_one (int thread_in_system) 
          {
          	         	
          	if (thread_in_system < nbod)
		{
			double x[3], v[3];
			_sys[thread_in_system].get(x[0],x[1],x[2],v[0],v[1],v[2]); // get current position and velocity
			
			BBox[6*thread_in_system] = min(X0[3*thread_in_system],x[0]); // xmin
			BBox[6*thread_in_system+1] = max(X0[3*thread_in_system],x[0]); // xmax
			
			BBox[6*thread_in_system+2] = min(X0[3*thread_in_system+1],x[1]); // ymin
			BBox[6*thread_in_system+3] = max(X0[3*thread_in_system+1],x[1]); // ymax
			
			BBox[6*thread_in_system+4] = min(X0[3*thread_in_system+2],x[2]); // zmin
			BBox[6*thread_in_system+5] = max(X0[3*thread_in_system+2],x[2]); // zmax
			
			// if velocity changes sign, do an interpolation
			
			for (int c = 0; c < 3; c++)
			{
			  if ( (V0[3*thread_in_system+c] < 0 && v[c] > 0 ) || 
				    (V0[3*thread_in_system+c] > 0 && v[c] < 0 ) )
			  {
				  double temp = V0[3*thread_in_system+c]*x[c] - v[c]*X0[3*thread_in_system+c]
							  - 0.5*timestep*V0[3*thread_in_system+c]*v[c]/(V0[3*thread_in_system+c]-v[c]);
				  BBox[6*thread_in_system + 2*c] = min(BBox[6*thread_in_system + 2*c],temp);
				  BBox[6*thread_in_system + 2*c+1] = max(BBox[6*thread_in_system + 2*c+1],temp);
				  
			  }
			}
			
			
			// Adjust values by the maximum close-encounter radius plus a fudge factor
			float temp = _sys[thread_in_system].attribute(0)*1.2; //Assume that rce is an attribute of each body
			for (int c = 0; c < 3; c++)
			{
			  BBox[6*thread_in_system + 2*c] -= temp;
			  BBox[6*thread_in_system + 2*c+1] += temp;
			}
		}
			
			
          }

        //! set the system state to disabled when three conditions are met
	GPUAPI void pass_two (int thread_in_system) 
      	{
      		if (thread_in_system < pair_count)
      		{
		      int i = nbod - 1 - thread_in_system / (nbod/2);
		      int j = thread_in_system % (nbod/2);
		      if (j >= i) 
		      {
			i = nbod - 1 - i - nbod%2 + 1;
			j = nbod - 1 - j - nbod%2 + 1;
		      }
		      // Check for close encounters between each pair of objects
	      
		      if ( BBox[6*i+1] >= BBox[6*j] && BBox[6*j+1] >= BBox[6*i] && 
			   BBox[6*i+3] >= BBox[6*j+2] && BBox[6*j+3] >= BBox[6*i+2] &&
			   BBox[6*i+5] >= BBox[6*j+4] && BBox[6*j+5] >= BBox[6*i+4] && 
			   _sys[i].mass() > 0 && _sys[j].mass() > 0)
		      {
			      float dx0 = X0[3*i] - X0[3*j];
			      float dy0 = X0[3*i+1] - X0[3*j+1];
			      float dz0 = X0[3*i+2] - X0[3*j+2];
			      float du0 = V0[3*i] - V0[3*j];
			      float dv0 = V0[3*i+1] - V0[3*j+1];
			      float dw0 = V0[3*i+2] - V0[3*j+2];
			      float d0t = (dx0*du0 + dy0*dv0 + dz0*dw0)*2.0; 
			      
			      float dx1 = _sys[i][0].pos() - _sys[j][0].pos();
			      float dy1 = _sys[i][1].pos() - _sys[j][1].pos();
			      float dz1 = _sys[i][2].pos() - _sys[j][2].pos();
			      float du1 = _sys[i][0].vel() - _sys[j][0].vel();
			      float dv1 = _sys[i][1].vel() - _sys[j][1].vel();
			      float dw1 = _sys[i][2].vel() - _sys[j][2].vel();
			      float d1t = (dx1*du1 + dy1*dv1 + dz1*dw1)*2.0;
			      
			      float d0 = dx0*dx0 + dy0*dy0 + dz0*dz0;
			      float d1 = dx1*dx1 + dy1*dy1 + dz1*dz1;
			      
			      float d2min, tmin, tau;
			      // Calculating the mininum distance between two trajectories of particle i and j
			      if (d0t > 0 && d1t < 0)
				      if (d0 < d1)
				      {
					      d2min = d0;
					      tmin = -timestep;
				      }
				      else
				      {
					      d2min = d1;
					      tmin = 0.0;
				      }
			      else // using Hermite interpolation to calculate trajectories and the possible intersections
			      {
				      float tmp = 6.0*(d0 - d1);
				      float a = tmp + 3.0*timestep*(d0t + d1t);
				      float b = tmp + 2.0*timestep*(d0t + 2.0*d1t);
				      float c = timestep * d1t;
				      
				      tmp = -0.5*(b + ((b>0)?1.0:(-1.0))*sqrt(max(b*b-4*a*c,0.0)));
				      if (tmp == 0.0) 
					tau = 0.0;
				      else
					tau = c/tmp;
				      tau = min(tau,0.0);
				      tau = max(tau,-1.0);
				      tmin = tau*timestep;
				      tmp = 1.0 + tau;
				      d2min = tau*tau*((3.0+2.0*tau)*d0 + tmp*timestep*d0t)
						      + tmp*tmp*((1.0-2.0*tau)*d1 + tau*timestep*d1t);
				      d2min = max(d2min, 0.0);
				      
			      }
			      
			      //---------------- End of minimum distance calculation --------------------------------
			      float d2ce = max(_sys[i].attribute(0), _sys[j].attribute(0)); // distant for closed encounter
			      float d2hit = max(_sys[i].attribute(1), _sys[j].attribute(1)); // distant for collision
			      d2ce = d2ce*d2ce;
			      d2hit = d2hit*d2hit;
			      float d2near = d2hit*4.0;
			      
			      if ((d2min <= d2ce && d0t <= 0 && d1t >= 0) || (d2min <= d2hit))
			      {
				/// do something with logging, even halting the system
			      }
			      
			      if (d2min < d2near)
			      {
				if (d2min < d2hit)
				{
				  /// Resolve collision
				  
				}
				/// Resolve closed encounter
			      }
			      
			      
		      }// endif: two bounding boxes have non-zero intersection
		} 
		/*
		if(is_condition_met() && is_deactivate_on() && (thread_in_system==0) )
		  {  _sys.set_disabled(); }
		return _sys.state();
		*/
	}
	
};

} // namespace monitor

} // namespace swarm
 
 
