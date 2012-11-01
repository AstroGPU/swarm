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

/*! \file mvs_cpu.hpp
 *   \brief Defines and implements \ref swarm::cpu::mvs_cpu class - the CPU 
 *          implementation of mixed variables symplectic propagator.
 *
 */

#ifndef H_MVS_CPU
#define H_MVS_CPU

#include "swarm/common.hpp"
#include "swarm/integrator.hpp"
#include "swarm/plugin.hpp"

//! Flag for using standard coordiates
#define  ASSUME_PROPAGATOR_USES_STD_COORDINATES 0

namespace swarm { namespace cpu {
/*! CPU implementation of mixed variables symplectic propagator: template<class Monitor>
 * \ingroup integrators
 *
 *   This is used as a reference implementation to
 *   test the GPU implementation of the mvs integrator
 *   
 *   This integrator can be used as an example of CPU integrator
 *
 *
 * \todo make Gravitation class a template parameter: template<class Monitor, class GravClass>
 */

template< class Monitor >
class mvs_cpu : public integrator {
	typedef integrator base;
	typedef Monitor monitor_t;
	typedef typename monitor_t::params mon_params_t;
	private:
	double _time_step;
	mon_params_t _mon_params;

  // included here so as to avoid namespace conflicts between CPU and OMP integrators
#include "../propagators/keplerian.hpp"

public:  //! Construct for class mvs_cpu
	mvs_cpu(const config& cfg): base(cfg),_time_step(0.001), _mon_params(cfg) {
		_time_step =  cfg.require("time_step", 0.0);
	}

	virtual void launch_integrator() {
		for(int i = 0; i < _ens.nsys(); i++){
			integrate_system(_ens[i]);
		}
	}

        //! Method for calculating inner product of two arrays
	inline static double inner_product(const double a[3],const double b[3]){
		return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
	}

        //! Method for calculating forces
	void calcForces(ensemble::SystemRef& sys, double acc[][3]){
		const int nbod = sys.nbod();

		// Clear acc 
		for(int b = 0; b < nbod; b++)	for(int c =0; c < 3; c++) 
		     acc[b][c] = 0.;

		// Loop through all pairs
		for(int i=0; i < nbod-1; i++) for(int j = i+1; j < nbod; j++) {

			double dx[3] = { sys[j][0].pos()-sys[i][0].pos(),
				sys[j][1].pos()-sys[i][1].pos(),
				sys[j][2].pos()-sys[i][2].pos()
			};
			double dv[3] = { sys[j][0].vel()-sys[i][0].vel(),
				sys[j][1].vel()-sys[i][1].vel(),
				sys[j][2].vel()-sys[i][2].vel()
			};

			// Calculated the magnitude
			double r2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2] * dx[2];
			double rinv = 1. / ( sqrt(r2) * r2 ) ;
			double rv =  inner_product(dx,dv) * 3. / r2;

			// Update acc for i
			const double scalar_i = +rinv*sys[j].mass();
			for(int c = 0; c < 3; c++) {
				acc[i][c] += dx[c]* scalar_i;
			}

			// Update acc for j
			const double scalar_j = -rinv*sys[i].mass();
			for(int c = 0; c < 3; c++) {
				acc[j][c] += dx[c]* scalar_j;
			}
		}
	}


  // MVS specific routines

	/// Standardized member name to call convert_helio_pos_bary_vel_to_std_coord 
	void convert_internal_to_std_coord(ensemble::SystemRef sys) 
	{ convert_helio_pos_bary_vel_to_std_coord(sys);	} 

  ///	Standardized member name to call convert_std_to_helio_pos_bary_vel_coord()
        GPUAPI void convert_std_to_internal_coord(ensemble::SystemRef sys) 
	{ convert_std_to_helio_pos_bary_vel_coord(sys); }


	/// Shift to  funky coordinate system (see A. Quillen's qymsym's tobary)
	void convert_std_to_helio_pos_bary_vel_coord(ensemble::SystemRef sys)  { 
	  const int nbod = sys.nbod();
	        double pc0;
		double sump = 0., sumv = 0., mtot = 0.;
		for(int c=0;c<3;++c)
		  {
		    pc0 = sys[0][c].pos();
		    // Find Center of mass and momentum
		    for(int j=0;j<nbod;++j) {
		      const double mj = sys[j].mass();
		      mtot += mj;
		      sump += mj*sys[j][c].pos();
		      sumv += mj*sys[j][c].vel();
		    }
		    sumv /= mtot;
		    int b = 0; // For sun
		    sys[b][c].vel() = sumv;
		    sys[b][c].pos() = sump/mtot;

		    for(b=1;b<nbod;++b) // For planets
		      {
			    sys[b][c].vel() -= sumv;
			    sys[b][c].pos() -= pc0;
		      }
		  }
	}

	/// Shift back from funky coordinate system (see A. Quillen's qymsym's frombary)
	void convert_helio_pos_bary_vel_to_std_coord (ensemble::SystemRef sys)  
	{ 
	  const int nbod = sys.nbod();
	  double sump = 0., sumv = 0., mtot;
	  double m0 = sys[0].mass();
	  double pc0, vc0;

	  for(int c=0;c<3;++c)
	    {
	      pc0 = sys[0][c].pos();
	      vc0 = sys[0][c].vel();
	      mtot = m0;
	      for(int j=1;j<nbod;++j)
		{
		  const double mj = sys[j].mass();
		  mtot += mj;
		  sump += mj*sys[j][c].pos();
		  sumv += mj*sys[j][c].vel();
		}
	      sump /= mtot;
	  
	      int b = 0; // For sun only
	      sys[b][c].pos() -= sump;
	      sys[b][c].vel() -= sumv/m0;

	      for(int b=1;b<nbod;++b) // For planets
		{
		  sys[b][c].pos() += pc0 - sump;
		  sys[b][c].vel() += vc0;
		}
	      
	    }
	}
	  

	/// Drift step for MVS integrator
  void drift_step(ensemble::SystemRef sys, const double hby2) 
	{
	  const double hby2m = hby2/sys[0].mass();
	  const int nbod = sys.nbod();
	  for(int c=0;c<3;++c)
	    {
	      int b = 0;
	      sys[b][c].pos() += hby2*sys[b][c].vel();
	      double mv = 0;
	      for(int j=1;j<nbod;++j)
		mv += sys[j].mass() * sys[j][c].vel();
	      for(b=1;b<nbod;++b)
		{	
		  sys[b][c].pos() += mv*hby2m;
		}
	    }
	}

        //! Integrating an ensemble
	void integrate_system(ensemble::SystemRef sys){
		const int nbod = sys.nbod();
		double acc[nbod][3];

		// Setting up Monitor
		monitor_t montest(_mon_params,sys,*_log) ;

		// begin init();
		const double sqrtGM = sqrt(sys[0].mass());
		convert_std_to_helio_pos_bary_vel_coord(sys);
		calcForces(sys,acc);
		// end init()

		for(int iter = 0 ; (iter < _max_iterations) && sys.is_active() ; iter ++ )
		  {

		// begin advance();
		    double hby2 = 0.5 * std::min( _destination_time - sys.time() ,  _time_step );

		// Step 1
		drift_step(sys,hby2);

		// Step 2: Kick Step
		for(int b=1;b<nbod;++b)
		  for(int c=0;c<3;++c)
		    sys[b][c].vel() += hby2 * acc[b][c];

		// __syncthreads();

		// 3: Kepler Drift Step (Keplerian orbit about sun/central body)
		for(int b=1;b<nbod;++b)
		  drift_kepler( sys[b][0].pos(),sys[b][1].pos(),sys[b][2].pos(),sys[b][0].vel(),sys[b][1].vel(),sys[b][2].vel(),sqrtGM, 2.0*hby2 );
		// __syncthreads();

		// TODO: check for close encounters here
		calcForces(sys,acc);

		// Step 4: Kick Step
		for(int b=1;b<nbod;++b)
		  for(int c=0;c<3;++c)
		    sys[b][c].vel() += hby2 * acc[b][c];
		// __syncthreads();

		// Step 5
		drift_step(sys,hby2);

		sys.time() += 2.0*hby2;

		// end advance

		const int thread_in_system_for_monitor = 0;
#if ASSUME_PROPAGATOR_USES_STD_COORDINATES
		montest( thread_in_system_for_monitor );
#else
		bool using_std_coord = false;
		bool needs_std_coord = montest.pass_one( thread_in_system_for_monitor );

		if(needs_std_coord) 
		  { 
		    convert_internal_to_std_coord(sys); 
		    using_std_coord = true; 
		  }
		
		//			__syncthreads();
		int new_state = montest.pass_two ( thread_in_system_for_monitor );

		if( montest.need_to_log_system() )
		  { log::system(*_log, sys); }
		
		//			__syncthreads();
		if(using_std_coord)
		  {
		    convert_std_to_internal_coord(sys);
		    using_std_coord = false;
		  }
#endif
			//			__syncthreads();			  

			if( sys.is_active() )
			  {
			    if( sys.time() >= _destination_time ) 
			      { sys.set_inactive();     }
			  }
			//			__syncthreads();

		  }

		// shutdown();
		convert_helio_pos_bary_vel_to_std_coord (sys);
	}
};



  } } // Close namespaces


#endif
