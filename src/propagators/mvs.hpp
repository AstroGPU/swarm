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

/*! \file mvs.hpp
 *   \brief Defines \ref swarm::gpu::bppt::MVSPropagator - the GPU implementation 
 *          of mixed variables symplectic propagator.
 *
 */

#include "swarm/common.hpp"
#include "swarm/swarmplugin.h"
#include "keplerian.hpp"

namespace swarm {

namespace gpu {
namespace bppt {

/*! Paramaters for MvsPropagator
 * \ingroup propagator_parameters
 *
 */
struct MVSPropagatorParams {
	double time_step;
        //! Constructor for MVSPropagatorParams
	MVSPropagatorParams(const config& cfg){
		time_step = cfg.require("time_step", 0.0);
	}
};

/*! GPU implementation of mixed variables symplectic propagator
 * \ingroup propagators
 *
 * \todo make Gravitation class a template parameter: template<class T, class GravClass>
 */
template<class T,class Gravitation>
struct MVSPropagator {
	typedef MVSPropagatorParams params;
	static const int nbod = T::n;

	params _params;


	//! Runtime variables
	ensemble::SystemRef& sys;
	Gravitation& calcForces;
	int b;
	int c;
	int ij;

	double sqrtGM;
	double max_timestep;

	double acc_bc;

        //! Constructor for MVSPropagator
	GPUAPI MVSPropagator(const params& p,ensemble::SystemRef& s,
			Gravitation& calc)
		:_params(p),sys(s),calcForces(calc){}

	__device__ bool is_in_body_component_grid()
        { return  ((b < nbod) && (c < 3)); }	

	__device__ bool is_in_body_component_grid_no_star()
        { return ( (b!=0) && (b < nbod) && (c < 3) ); }	

	__device__ bool is_first_thread_in_system()
        { return (thread_in_system()==0); }	

	static GENERIC int thread_per_system(){
		return nbod * 3;
	}

	static GENERIC int shmem_per_system() {
		 return 0;
	}



	/// Shift into funky coordinate system (see A. Quillen's qymsym's tobary)
	/// Shift back and forth is tested and it is indeed symmetric 
	/// Initialization tasks executed before entering loop
        /// Cache sqrtGM, shift coord system, cache acceleration data for this thread's body and component
	GPUAPI void init()  { 
		sqrtGM = sqrt(sys[0].mass());
		convert_std_to_helio_pos_bary_vel_coord();
		__syncthreads();
		acc_bc = calcForces.acc_planets(ij,b,c);
                }

	/// Before exiting, convert back to standard cartesian coordinate system
	GPUAPI void shutdown() { 
	convert_helio_pos_bary_vel_to_std_coord ();
	}

	/// Shift to  funky coordinate system (see A. Quillen's qymsym's tobary)
	GPUAPI void convert_std_to_helio_pos_bary_vel_coord()  
	{ convert_std_to_helio_pos_bary_vel_coord_without_shared();  }

	/// Shift to  funky coordinate system (see A. Quillen's qymsym's tobary)
	/// At least when there's not logging that needs frequent shifts, this is too small to bother with shared memory
	GPUAPI void convert_std_to_helio_pos_bary_vel_coord_with_shared()  { 
		double sump = 0., sumv = 0., mtot = 0.;
		if( is_in_body_component_grid() )
		{

		if( b==0 )
		{
			calcForces.shared[1][c].acc() = sys[0][c].pos();
			//! Find Center of mass and momentum
			for(int j=0;j<nbod;++j) {
				const double mj = sys[j].mass();
				mtot += mj;
				sump += mj*sys[j][c].pos();
				sumv += mj*sys[j][c].vel();
			}
		sumv /= mtot;
		// There should be a better way to access shared memory
		// What if we need to change coordinates
		calcForces.shared[0][c].acc() = sumv;
		}

		}
		__syncthreads();

		if( is_in_body_component_grid() )
		{
			if(b==0) // For sun
			{
			sys[b][c].vel() = sumv;
			sys[b][c].pos() = sump/mtot;
			}
			else     // For planets
			{
			sys[b][c].vel() -= calcForces.shared[0][c].acc(); // really sumv from shared;
			sys[b][c].pos() -= calcForces.shared[1][c].acc(); // really pc0 = original sys[0][c].pos() from shared
			}
		}

	}

	
	/// Shift to  funky coordinate system (see A. Quillen's qymsym's tobary)
	GPUAPI void convert_std_to_helio_pos_bary_vel_coord_without_shared()  { 
	        double pc0;
		double sump = 0., sumv = 0., mtot = 0.;
		if( is_in_body_component_grid() )
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
		}
		__syncthreads();

		if( is_in_body_component_grid() )
		{
			if(b==0) // For sun
			{
			sys[b][c].vel() = sumv;
			sys[b][c].pos() = sump/mtot;
			}
			else     // For planets
			{
			sys[b][c].vel() -= sumv;
			sys[b][c].pos() -= pc0;
			}
		}

	}



	/// Shift back from funky coordinate system (see A. Quillen's qymsym's frombary)
	GPUAPI void convert_helio_pos_bary_vel_to_std_coord ()  
	{ 
	  double sump = 0., sumv = 0., mtot;
	  double m0, pc0, vc0;
		if ( is_in_body_component_grid() ) 
		{
		  m0 = sys[0].mass();
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
		}
		__syncthreads();

		if ( is_in_body_component_grid() ) 
		   {
		   if(b==0) // For sun only
			{
			sys[b][c].pos() -= sump;
			sys[b][c].vel() -= sumv/m0;
			}
		   else
			{
			sys[b][c].pos() += pc0 - sump;
			sys[b][c].vel() += vc0;
			}
		}

	}

	/// Standardized member name to call convert_helio_pos_bary_vel_to_std_coord 
	GPUAPI void convert_internal_to_std_coord() 
	{ convert_helio_pos_bary_vel_to_std_coord ();	} 

	/// Standardized member name to call convert_std_to_helio_pos_bary_vel_coord_without_shared()
        GPUAPI void convert_std_to_internal_coord() 
	{ convert_std_to_helio_pos_bary_vel_coord_without_shared(); }


	/// Drift step for MVS integrator
	GPUAPI void drift_step(const double hby2) 
	{
	      if(b==0)
	        {
		sys[b][c].pos() += hby2*sys[b][c].vel();
		}
	      else
	        {
	      	double mv = 0;
	      	for(int j=1;j<nbod;++j)
	      		mv += sys[j].mass() * sys[j][c].vel();

		sys[b][c].pos() += mv*hby2/sys[0].mass();
		}
	}


	/// Advance system by one time unit
	GPUAPI void advance()
	{
		double hby2 = 0.5 * min( max_timestep ,  _params.time_step );

			// Step 1
			if ( is_in_body_component_grid() ) 
			   drift_step(hby2);

			// Step 2: Kick Step
			if( is_in_body_component_grid_no_star() ) 
			   sys[b][c].vel() += hby2 * acc_bc;

			__syncthreads();

			// 3: Kepler Drift Step (Keplerian orbit about sun/central body)
			if( (ij>0) && (ij<nbod)  ) 
			    drift_kepler( sys[ij][0].pos(),sys[ij][1].pos(),sys[ij][2].pos(),sys[ij][0].vel(),sys[ij][1].vel(),sys[ij][2].vel(),sqrtGM, 2.0*hby2 );
			__syncthreads();

			// TODO: check for close encounters here
			acc_bc = calcForces.acc_planets(ij,b,c);

			// Step 4: Kick Step
			if( is_in_body_component_grid_no_star() ) 
			   sys[b][c].vel() += hby2 * acc_bc;
			__syncthreads();

			// Step 5
			if ( is_in_body_component_grid() ) 
			  drift_step(hby2);

		if( is_first_thread_in_system() ) 
			sys.time() += 2.0*hby2;
	}
};




}
}
}

