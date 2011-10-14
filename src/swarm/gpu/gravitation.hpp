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
#include "../types/coalescedstructarray.hpp"


namespace swarm {
namespace gpu {
namespace bppt {

/**
 *  Unit type of the acceleration and jerk pairs shared array.
 *   
 *  This for each pair we keep an acc and a jerk. The values are
 *  not final values, they are intermediate values calculated by
 *  calc_pair and should be accumulated using the correct algorithm
 *  to produce acceleration and jerk. WARPSIZE can be 1. Usually it is
 *  set to 16 for optimizing coalesced reads from memory.
 *
 */
template<int W>
struct GravitationScalars {
	static const int WARPSIZE = W;
	typedef double scalar_t; 

	double _acc[WARPSIZE];
	double _jerk[WARPSIZE];

	// Accessors
	GENERIC double& acc() { return _acc[0];  }
	GENERIC double& jerk() { return _jerk[0];  }
};


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
template<int nbod, int WARPSIZE = SHMEM_WARPSIZE>
class Gravitation {
	public:
	const static int pair_count = (nbod*(nbod-1))/2;

	typedef GravitationScalars<WARPSIZE> shared_data [pair_count][3];

	private:
	ensemble::SystemRef& sys;
	shared_data &shared;

	/**
	 * Helper function for calculating inner product
	 */
	GENERIC static double inner_product(const double a[3],const double b[3]){
		return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
	}

	/// Helper function to convert an integer from 1..n*(n-1)/2 to a pair (first,second), this function returns the first element.
	GENERIC static int first ( int ij ){
		int i = nbod - 1 - ij / (nbod/2);
		int j = ij % (nbod/2);
		if (j < i) 
			return i;
		else 
			return nbod - 1 - i - nbod%2 + 1;
	}

	/// Helper function to convert an integer from 1..n*(n-1)/2 to a pair (first,second), this function returns the second element.
	GENERIC static int second ( int ij ){
		int i = nbod - 1 - ij / (nbod/2);
		int j = ij % (nbod/2);
		if (j < i) 
			return j;
		else 
			return nbod - 1 - j - nbod%2;
	}

	public:

	/*
	 * Create a function object for computing gravitational force 
	 * on planets in a system using a shared memory area.
	 *
	 * @arg sys   Reference to system that this algorithm operates on
	 * @arg shared Reference to an array of appropriate size
	 * allocated on shared memory to hold intermediat results.
	 *
	 */
	GENERIC Gravitation(ensemble::SystemRef& sys,shared_data &shared):sys(sys),shared(shared){	}

	/**
	 *  Step one of the algorithm. All pairs run in parallel. This
	 *  function calculates intermediate results for a pair and
	 *  stores it into shared array.
	 *
	 *  @arg ij Integer number refering to a pair
	 */
	GENERIC void calc_pair(int ij)const{
		int i = first( ij );
		int j = second( ij );
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

	/*
	 * calculate accleration for a planet ignoring the
	 * impact of the body 0 (star).
	 *  @b  planet number
	 *  @c  coordinate number x:0,y:1,z:2
	 */
	GENERIC double sum_acc_planets(int b,int c)const{
		double total = 0;

#pragma unroll
		for(int d = 0; d < pair_count; d++){
			int x = first(d), y= second(d);

			if(x == b){
				if(y != 0)
					total += shared[d][c].acc() * sys[y].mass();
			}else if(y == b){
				if(x != 0)
					total -= shared[d][c].acc() * sys[x].mass();
			}
		}

		return total;
	}

	/*  
	 *  Find the acceleration for a planet.
	 *
	 *  @b  planet number
	 *  @c  coordinate number x:0,y:1,z:2
	 */
	GENERIC double sum_acc(int b,int c)const{
		double total = 0;
		double from_sun = 0;

		/// Find the contribution from/to Sun first
#pragma unroll
		for(int d = 0; d < pair_count; d++){
			int x = first(d), y= second(d);

			if(x == b){
				if(y == 0)
					from_sun += shared[d][c].acc() * sys[y].mass();
				else
					total += shared[d][c].acc() * sys[y].mass();
			}else if(y == b){
				if(x == 0)
					from_sun -= shared[d][c].acc() * sys[x].mass();
				else
					total -= shared[d][c].acc() * sys[x].mass();
			}
		}

		return from_sun + total;
	}

	/*  
	 *  Find the acceleration and jerk for a planet.
	 *
	 *  @b  planet number
	 *  @c  coordinate number x:0,y:1,z:2
	 *  @acc  reference to output variable for acceleration
	 *  @jerk reference to output variable for jerk
	 */
	GENERIC void sum(int b,int c,double& acc, double & jerk)const{
		// Total acceleration from other planets
		double total_a = 0.0;
		// Total acceleration from body 0 (sun or star)
		double from_sun_a = 0.0;
		// Total jerk from other planets
		double total_j = 0.0;
		// Total jerk from body 0 (sun or star)
		double from_sun_j = 0.0;

		/// Find the contribution from/to Sun first
#pragma unroll
		for(int d = 0; d < pair_count; d++){
			int x = first(d), y= second(d);

			if(x == b){
				if(y == 0) {
					from_sun_a += shared[d][c].acc() * sys[y].mass();
					from_sun_j += shared[d][c].jerk() * sys[y].mass();
				} else {
					total_a += shared[d][c].acc() * sys[y].mass();
					total_j += shared[d][c].jerk() * sys[y].mass();
				}
			}else if(y == b){
				if(x == 0) {
					from_sun_a -= shared[d][c].acc() * sys[x].mass();
					from_sun_j -= shared[d][c].jerk() * sys[x].mass();
				} else {
					total_a -= shared[d][c].acc() * sys[x].mass();
					total_j -= shared[d][c].jerk() * sys[x].mass();
				}
			}
		}

		acc = from_sun_a + total_a;
		jerk = from_sun_j + total_j;
	}


	/*
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
	GPUAPI void operator() (int ij,int b,int c,double& pos,double& vel,double& acc,double& jerk)const{
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

	/*
	 * Different version of acceleration calculation used for 
	 * MVS integrator. The impact of body 0(sun or star) is 
	 * ignored because in the integrator it is calculated using
	 * keplerian motion.
	 * This is tightly coupled with the
	 * BPPT integrators. ij, b and c are calculated from thread id.
	 *
	 * @ij The pair number for this tread.
	 * @b  The planet number for this thread.
	 * @c  coordinate number x:0,y:1,z:2
	 * @pos position for this planet's coordinate
	 * @vel velecotiy for this planet's coordinate
	 *
	 */
	GPUAPI double acc_planets (int ij,int b,int c)const{
		if(ij < pair_count)
			calc_pair(ij);
		__syncthreads();
		if(b < nbod && c < 3){
			return sum_acc_planets(b,c);
		}else
			return 0;
	}

	/*
	 * Run the complete algorithm for computing acceleration only 
	 * on all bodies. This is tightly coupled with the
	 * BPPT integrators. ij, b and c are calculated from thread id.
	 *
	 * @ij The pair number for this tread.
	 * @b  The planet number for this thread.
	 * @c  coordinate number x:0,y:1,z:2
	 * @pos position for this planet's coordinate
	 * @vel velecotiy for this planet's coordinate
	 *
	 */
	GPUAPI double acc (int ij,int b,int c,double& pos,double& vel)const{
		// Write positions to shared (global) memory
		if(b < nbod && c < 3)
			sys[b][c].pos() = pos , sys[b][c].vel() = vel;
		__syncthreads();
		if(ij < pair_count)
			calc_pair(ij);
		__syncthreads();
		if(b < nbod && c < 3){
			return sum_acc(b,c);
		}else
			return 0;
	}


};

}
}
}
