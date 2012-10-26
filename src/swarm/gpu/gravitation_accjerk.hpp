/*************************************************************************
 * Copyright (C) 2010 by Saleh Dindar and the Swarm-NG Development Team  *
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

/*! \file gravitation_accjerk.hpp
 *   \brief Defines and implements class \ref swarm::gpu::bppt::GravitationAccJerk
 *          that implements the functions to calculate accerleration and jerk 
 *          of the gravitation in parallel.
 *          
 */

#pragma once

#include "gravitation_common.hpp"

namespace swarm { namespace gpu { namespace bppt {

/** 
 * templatized Class working as a function object to 
 * calculate acceleration and jerk in parallel.
 *
 * shared_data member is to be initialized with a shared memory pointer.
 * Although it does not have to be a shared memory pointer.
 *
 * It operates in two steps:
 *
 * Step 0: Write positions and velocities to global memory which is cached. 
 *
 * Step 1: Calculate distances between pairs using \ref calc_pair
 * you should supply ij that is between 0 and n*(n-1)/2. It calculates the
 * inverse of distance squared between each pair of bodies. The
 * intermediate values for calculating acc/jerk are stored to the
 * shared memory.
 *
 * Step 2: Calculate forces that operate on each body it is done
 * in one of the functions: sum, sum_acc, sum_acc_planets
 * This function should be called per body per component. It uses the shared
 * data and calculates the acc/jerk for each body.
 *
 */
template<class T>
class GravitationAccJerk {
	public:
	const static int nbod = T::n;
	const static int pair_count = (nbod*(nbod-1))/2;
	const static int CHUNK_SIZE = SHMEM_CHUNK_SIZE;

	typedef GravitationAccJerkScalars<CHUNK_SIZE> shared_data [pair_count][3];

	private:
	public: // hack to make this public to test whether it's worth using shared memory for some small steps
	ensemble::SystemRef& sys;
	shared_data &shared;


	public:

	/**
	 * Create a function object for computing gravitational force 
	 * on planets in a system using a shared memory area.
	 *
	 * @arg sys   Reference to system that this algorithm operates on
	 * @arg shared Reference to an array of appropriate size
	 * allocated on shared memory to hold intermediat results.
	 *
	 */
	GENERIC GravitationAccJerk(ensemble::SystemRef& sys,shared_data &shared):sys(sys),shared(shared){	}

	private:

	/**
	 *  Step one of the algorithm. All pairs run in parallel. This
	 *  function calculates intermediate results for a pair and
	 *  stores it into shared array.
	 *
	 *  @arg ij Integer number refering to a pair
	 */
	GENERIC void calc_pair(int ij)const{
		int i = first<nbod>( ij );
		int j = second<nbod>( ij );
		if(i != j){

			// Relative vector from planet i to planet j
			double dx[3] =  { sys[j][0].pos()- sys[i][0].pos(),sys[j][1].pos()- sys[i][1].pos(), sys[j][2].pos()- sys[i][2].pos() };
			// Relative velocity between the planets
			double dv[3] =  { sys[j][0].vel()- sys[i][0].vel(),sys[j][1].vel()- sys[i][1].vel(), sys[j][2].vel()- sys[i][2].vel() };

			// Distance between the planets
			double r2 =  dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2];

			double jerk_mag =  inner_product(dx,dv) * 3.0 / r2;
			double acc_mag =  rsqrt(r2) / r2;

#pragma unroll
			for(int c = 0; c < 3; c ++)
			{
				shared[ij][c].acc() = dx[c]* acc_mag;
				shared[ij][c].jerk() = (dv[c] - dx[c] * jerk_mag ) * acc_mag;
			}
		}

	}


	public:

	/**  
	 *  Find the acceleration and jerk for a planet.
	 *
	 *  @b  planet number
	 *  @c  coordinate number x:0,y:1,z:2
	 *  @acc  reference to output variable for acceleration
	 *  @jerk reference to output variable for jerk
	 */
	GENERIC void sum(int b,int c,double& acc, double & jerk)const{
		// Total acceleration from other planets
		double acc_from_planets = 0.0;
		// Total acceleration from body 0 (sun or star)
		double acc_from_sun = 0.0;
		// Total jerk from other planets
		double jerk_from_planets = 0.0;
		// Total jerk from body 0 (sun or star)
		double jerk_from_sun = 0.0;

		/// Find the contribution from/to Sun first
#pragma unroll
		for(int d = 0; d < pair_count; d++){
			int x = first<nbod>(d), y= second<nbod>(d);

			if(x == b){
				if(y == 0) {
					acc_from_sun += shared[d][c].acc() * sys[y].mass();
					jerk_from_sun += shared[d][c].jerk() * sys[y].mass();
				} else {
					acc_from_planets += shared[d][c].acc() * sys[y].mass();
					jerk_from_planets += shared[d][c].jerk() * sys[y].mass();
				}
			}else if(y == b){
				if(x == 0) {
					acc_from_sun -= shared[d][c].acc() * sys[x].mass();
					jerk_from_sun -= shared[d][c].jerk() * sys[x].mass();
				} else {
					acc_from_planets -= shared[d][c].acc() * sys[x].mass();
					jerk_from_planets -= shared[d][c].jerk() * sys[x].mass();
				}
			}
		}

		acc = acc_from_sun + acc_from_planets;
		jerk = jerk_from_sun + jerk_from_planets;
	}

	public:

	/**
	 * Run the complete algorithm for computing acceleration and
	 * jerk on all bodies. This is tightly coupled with the
	 * BPPT integrators. ij, b and c are calculated from thread id.
	 *
	 * If you need to calculate only acceleration use \ref acc function
	 * instead.
	 *
	 * @ij The pair number for this tread.
	 * @b  The planet number for this thread.
	 * @c  coordinate number x:0,y:1,z:2
	 * @pos position for this planet's coordinate
	 * @vel velecotiy for this planet's coordinate
	 * @acc output variable to hold acceleration
	 * @jerk output variable to hold jerk.
	 *
	 */
	GPUAPI void operator() (int ij,int b,int c,double& pos,double& vel,double& acc,double& jerk) const{
		// Write positions to shared (global) memory
		if(b < nbod && c < 3)
			sys[b][c].pos() = pos , sys[b][c].vel() = vel;
		__syncthreads();
		if(ij < pair_count)
			calc_pair(ij);
		__syncthreads();
		if(b < nbod && c < 3){
			sum(b,c,acc,jerk);
		}
	}
				
	static GENERIC int thread_per_system(){
		return (nbod*3>(nbod-1)*nbod/2) ? nbod*3 : (nbod-1)*nbod/2;
	}

	static GENERIC int shmem_per_system() {
		 return sizeof(shared_data)/CHUNK_SIZE;
	}


				
};

} } } // end of namespace bppt :: gpu :: swarm
