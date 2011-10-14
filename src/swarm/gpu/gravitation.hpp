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
#ifndef H_GRAVITATION
#define H_GRAVITATION
// TODO: Do we actually use this?
#include "../types/coalescedstructarray.hpp"
#include "bppt.hpp"

namespace swarm {
namespace gpu {
namespace bppt {

/**
 *  Unit type of the acceleration pairs shared array.
 *   
 *  This for each pair we keep an acc. The values are
 *  not final values, they are intermediate values calculated by
 *  calc_pair and should be accumulated using the correct algorithm
 *  to produce acceleration. WARPSIZE can be 1. Usually it is
 *  set to 16 for optimizing coalesced reads from memory.
 *
 */
template<int W>
struct GravitationAccOnlyScalars {
	static const int WARPSIZE = W;
	typedef double scalar_t; 

	double _acc[WARPSIZE];

	// Accessors
	GENERIC double& acc() { return _acc[0];  }
};


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


  // Default warpsize suggestion for large systems
  template<int nbod, int gpu_shared_mem_size>  struct SelectWarpSizeGravitation
  {    static const int suggest = 1;  };

  // Default warpsize suggestion for small systems (assuming 16k cache)
  template<>  struct SelectWarpSizeGravitation<2,16384>
  {    static const int suggest = 16;  };
  template<>  struct SelectWarpSizeGravitation<3,16384>
  {    static const int suggest = 16;  };
  template<>  struct SelectWarpSizeGravitation<4,16384>
  {    static const int suggest = 16;  };
  template<>  struct SelectWarpSizeGravitation<5,16384>
  {    static const int suggest = 16;  };
  template<>  struct SelectWarpSizeGravitation<6,16384>
  {    static const int suggest = 16;  };

  // Default warpsize suggestion for small systems (assuming 48k cache)
  template<>  struct SelectWarpSizeGravitation<2,49152>
  {    static const int suggest = 16;  };
  template<>  struct SelectWarpSizeGravitation<3,49152>
  {    static const int suggest = 16;  };
  template<>  struct SelectWarpSizeGravitation<4,49152>
  {    static const int suggest = 16;  };
  template<>  struct SelectWarpSizeGravitation<5,49152>
  {    static const int suggest = 16;  };
  template<>  struct SelectWarpSizeGravitation<6,49152>
  {    static const int suggest = 16;  };
  template<>  struct SelectWarpSizeGravitation<7,49152>
  {    static const int suggest = 16;  };
  template<>  struct SelectWarpSizeGravitation<8,49152>
  {    static const int suggest = 16;  };
  template<>  struct SelectWarpSizeGravitation<9,49152>
  {    static const int suggest = 16;  };
  template<>  struct SelectWarpSizeGravitation<10,49152>
  {    static const int suggest = 16;  };
  template<>  struct SelectWarpSizeGravitation<11,49152>
  {    static const int suggest = 16;  };

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
//template<int nbod, int WARPSIZE = SHMEM_WARPSIZE>
  template<int nbod, int WARPSIZE = SelectWarpSizeGravitation<nbod,MIN_SHMEM_SIZE>::suggest >
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
        __device__ Gravitation(ensemble::SystemRef& sys, const int  sysid_in_block):sys(sys),
							 shared(*( (shared_data*) system_shared_data_pointer( sysid_in_block)))
				{	  }

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

  // TODO: Remove once allow propagators to use GravitationAccOnly
	/*
	 * calculate accleration for a planet ignoring the
	 * impact of the body 0 (star).
	 *  @b  planet number
	 *  @c  coordinate number x:0,y:1,z:2
	 */
	GENERIC double sum_acc_planets(int b,int c)const{
		double acc_sum = 0;

#pragma unroll
		for(int d = 0; d < pair_count; d++){
			int x = first(d), y= second(d);

			if(x == b){
				if(y != 0)
					acc_sum += shared[d][c].acc() * sys[y].mass();
			}else if(y == b){
				if(x != 0)
					acc_sum -= shared[d][c].acc() * sys[x].mass();
			}
		}

		return acc_sum;
	}

  // TODO: Remove once allow propagators to use GravitationAccOnly
	/*  
	 *  Find the acceleration for a planet.
	 *
	 *  @b  planet number
	 *  @c  coordinate number x:0,y:1,z:2
	 */
	GENERIC double sum_acc(int b,int c)const{
		double acc_from_planets = 0;
		double acc_from_sun = 0;

		/// Find the contribution from/to Sun first
#pragma unroll
		for(int d = 0; d < pair_count; d++){
			int x = first(d), y= second(d);

			if(x == b){
				if(y == 0)
					acc_from_sun += shared[d][c].acc() * sys[y].mass();
				else
					acc_from_planets += shared[d][c].acc() * sys[y].mass();
			}else if(y == b){
				if(x == 0)
					acc_from_sun -= shared[d][c].acc() * sys[x].mass();
				else
					acc_from_planets -= shared[d][c].acc() * sys[x].mass();
			}
		}

		return acc_from_sun + acc_from_planets;
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
			int x = first(d), y= second(d);

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
				
  // TODO: Remove once allow propagators to use GravitationAccOnly
	__device__ double acc_planets (int ij,int b,int c)const{
		if(ij < pair_count)
			calc_pair(ij);
		__syncthreads();
		if(b < nbod && c < 3){
			return sum_acc_planets(b,c);
		}else
			return 0;
	}
				
	__device__ double acc (int ij,int b,int c,double& pos,double& vel)const{
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


	static GENERIC int shmem_per_system() {
		const int pair_count = nbod * (nbod - 1) / 2;
		//		return pair_count * 3  * 2 * sizeof(double);
		// TODO: Test
		 return pair_count * 3  * sizeof(GravitationScalars<WARPSIZE>)/WARPSIZE;
	}

	static __device__ void * system_shared_data_pointer(const int sysid_in_block) {
		extern __shared__ char shared_mem[];
		int b = sysid_in_block / WARPSIZE ;
		int i = sysid_in_block % WARPSIZE ;
		int idx = i * sizeof(double) 
			+ b * WARPSIZE 
			* shmem_per_system();
		return &shared_mem[idx];
	}

	static __device__ void * unused_shared_data_pointer(const int system_per_block) {
		extern __shared__ char shared_mem[];
		int idx = system_per_block * shmem_per_system();
		return &shared_mem[idx];
	}

				
};


  // Default warpsize suggestion for large systems
  template<int nbod, int gpu_shared_mem_size>  struct SelectWarpSizeGravitationAccOnly
  {    static const int suggest = 1;  };

  // Default warpsize suggestion for small systems (assuming 16k cache)
  template<>  struct SelectWarpSizeGravitationAccOnly<2,16384>
  {    static const int suggest = 16;  };
  template<>  struct SelectWarpSizeGravitationAccOnly<3,16384>
  {    static const int suggest = 16;  };
  template<>  struct SelectWarpSizeGravitationAccOnly<4,16384>
  {    static const int suggest = 16;  };
  template<>  struct SelectWarpSizeGravitationAccOnly<5,16384>
  {    static const int suggest = 16;  };
  template<>  struct SelectWarpSizeGravitationAccOnly<6,16384>
  {    static const int suggest = 16;  };
  template<>  struct SelectWarpSizeGravitationAccOnly<7,16384>
  {    static const int suggest = 16;  };
  template<>  struct SelectWarpSizeGravitationAccOnly<8,16384>
  {    static const int suggest = 16;  };
  template<>  struct SelectWarpSizeGravitationAccOnly<9,16384>
  {    static const int suggest = 16;  };

  // Default warpsize suggestion for small systems (assuming 48k cache)
  template<>  struct SelectWarpSizeGravitationAccOnly<2,49152>
  {    static const int suggest = 16;  };
  template<>  struct SelectWarpSizeGravitationAccOnly<3,49152>
  {    static const int suggest = 16;  };
  template<>  struct SelectWarpSizeGravitationAccOnly<4,49152>
  {    static const int suggest = 16;  };
  template<>  struct SelectWarpSizeGravitationAccOnly<5,49152>
  {    static const int suggest = 16;  };
  template<>  struct SelectWarpSizeGravitationAccOnly<6,49152>
  {    static const int suggest = 16;  };
  template<>  struct SelectWarpSizeGravitationAccOnly<7,49152>
  {    static const int suggest = 16;  };
  template<>  struct SelectWarpSizeGravitationAccOnly<8,49152>
  {    static const int suggest = 16;  };
  template<>  struct SelectWarpSizeGravitationAccOnly<9,49152>
  {    static const int suggest = 16;  };
  template<>  struct SelectWarpSizeGravitationAccOnly<10,49152>
  {    static const int suggest = 16;  };
  template<>  struct SelectWarpSizeGravitationAccOnly<11,49152>
  {    static const int suggest = 16;  };
  template<>  struct SelectWarpSizeGravitationAccOnly<12,49152>
  {    static const int suggest = 16;  };
  template<>  struct SelectWarpSizeGravitationAccOnly<13,49152>
  {    static const int suggest = 16;  };
  template<>  struct SelectWarpSizeGravitationAccOnly<14,49152>
  {    static const int suggest = 16;  };
  template<>  struct SelectWarpSizeGravitationAccOnly<15,49152>
  {    static const int suggest = 16;  };
  template<>  struct SelectWarpSizeGravitationAccOnly<16,49152>
  {    static const int suggest = 16;  };

/*! 
 * templatized Class to calculate acceleration and jerk in parallel
 *
 * Similar to Gravitation, but computes only terms for acceleration and not for jerk
 * To be used with integration aglorithms that don't make use of the jerk, e.g., Runge-Kutta
 */

  template<int nbod, int WARPSIZE = SelectWarpSizeGravitationAccOnly<nbod,MIN_SHMEM_SIZE>::suggest >
class GravitationAccOnly {
	public:
	const static int pair_count = (nbod*(nbod-1))/2;

	typedef GravitationAccOnlyScalars<WARPSIZE> shared_data [pair_count][3];

	private:
	ensemble::SystemRef& sys;
	shared_data &shared;

	inline __device__ static double inner_product(const double a[3],const double b[3]){
		return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
	}

	__device__ static int first ( int ij ){
		int i = nbod - 1 - ij / (nbod/2);
		int j = ij % (nbod/2);
		if (j < i) 
			return i;
		else 
			return nbod - 1 - i - nbod%2 + 1;
	}

	__device__ static int second ( int ij ){
		int i = nbod - 1 - ij / (nbod/2);
		int j = ij % (nbod/2);
		if (j < i) 
			return j;
		else 
			return nbod - 1 - j - nbod%2;
	}

	public:

	__device__ GravitationAccOnly(ensemble::SystemRef& sys,shared_data &shared):sys(sys),shared(shared){	}
  // TODO: Need to test this version
  __device__ GravitationAccOnly(ensemble::SystemRef& sys, int sysid_in_block):sys(sys),
							  shared(*( (shared_data*) system_shared_data_pointer(sysid_in_block))) {  }


	__device__ void calc_pair(int ij)const{
		int i = first( ij );
		int j = second( ij );
		if(i != j){

			double dx[3] =  { sys[j][0].pos()- sys[i][0].pos(),sys[j][1].pos()- sys[i][1].pos(), sys[j][2].pos()- sys[i][2].pos() };
			double dv[3] =  { sys[j][0].vel()- sys[i][0].vel(),sys[j][1].vel()- sys[i][1].vel(), sys[j][2].vel()- sys[i][2].vel() };
			double r2 =  dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2];
			double acc_mag =  rsqrt(r2) / r2;

#pragma unroll
			for(int c = 0; c < 3; c ++)
			{
				shared[ij][c].acc() = dx[c]* acc_mag;
			}
		}

	}

	__device__ double sum_acc_planets(int b,int c)const{
		double acc_sum = 0;

		/// Find the contribution from/to Sun first
#pragma unroll
		for(int d = 0; d < pair_count; d++){
			int x = first(d), y= second(d);

			if(x == b){
				if(y != 0)
					acc_sum += shared[d][c].acc() * sys[y].mass();
			}else if(y == b){
				if(x != 0)
					acc_sum -= shared[d][c].acc() * sys[x].mass();
			}
		}

		return acc_sum;
	}

	__device__ double sum_acc(int b,int c)const{
		double acc_from_planets = 0;
		double acc_from_sun = 0;

		/// Find the contribution from/to Sun first
#pragma unroll
		for(int d = 0; d < pair_count; d++){
			int x = first(d), y= second(d);

			if(x == b){
				if(y == 0)
					acc_from_sun += shared[d][c].acc() * sys[y].mass();
				else
					acc_from_planets += shared[d][c].acc() * sys[y].mass();
			}else if(y == b){
				if(x == 0)
					acc_from_sun -= shared[d][c].acc() * sys[x].mass();
				else
					acc_from_planets -= shared[d][c].acc() * sys[x].mass();
			}
		}

		return acc_from_sun + acc_from_planets;
	}


	__device__ void operator() (int ij,int b,int c,double& pos,double& vel,double& acc)const{
		// Write positions to shared (global) memory
		if(b < nbod && c < 3)
			sys[b][c].pos() = pos , sys[b][c].vel() = vel;
		__syncthreads();
		if(ij < pair_count)
			calc_pair(ij);
		__syncthreads();
		if(b < nbod && c < 3){
			acc = sum_acc(b,c);
		}
	}

	  // WARNING: Should this function copy pos and vel to global first?
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


	static GENERIC int shmem_per_system() {
		const int pair_count = nbod * (nbod - 1) / 2;
		return pair_count * 3  * sizeof(double);
		// TODO: Test
		// return pair_count * 3  * sizeof(GravitationScalars<WARPSIZE>)/WARPSIZE;
	}

	static __device__ void * system_shared_data_pointer(const int sysid_in_block) {
		extern __shared__ char shared_mem[];
		int b = sysid_in_block / WARPSIZE ;
		int i = sysid_in_block % WARPSIZE ;
		int idx = i * sizeof(double) 
			+ b * WARPSIZE 
			* shmem_per_system();
		return &shared_mem[idx];
	}

	static __device__ void * unused_shared_data_pointer(const int system_per_block) {
		extern __shared__ char shared_mem[];
		int idx = system_per_block * shmem_per_system();
		return &shared_mem[idx];
	}


}
;
}
}
}

#endif
