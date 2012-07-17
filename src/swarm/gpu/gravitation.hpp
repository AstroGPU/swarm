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
#pragma once

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
 *  to produce acceleration. CHUNK_SIZE can be 1. Usually it is
 *  set to 16 for optimizing coalesced reads from memory.
 *
 */
template<int W>
struct GravitationAccOnlyScalars {
	static const int CHUNK_SIZE = W;
	typedef double scalar_t; 

	double _acc[CHUNK_SIZE];

	// Accessors
	GENERIC double& acc() { return _acc[0];  }
};


/**
 *  Unit type of the acceleration and jerk pairs shared array.
 *   
 *  This for each pair we keep an acc and a jerk. The values are
 *  not final values, they are intermediate values calculated by
 *  calc_pair and should be accumulated using the correct algorithm
 *  to produce acceleration and jerk. CHUNK_SIZE can be 1. Usually it is
 *  set to 16 for optimizing coalesced reads from memory.
 *
 */
template<int W>
struct GravitationScalars {
	static const int CHUNK_SIZE = W;
	typedef double scalar_t; 

	double _acc[CHUNK_SIZE];
	double _jerk[CHUNK_SIZE];

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
template<int nbod, int CHUNK_SIZE = SHMEM_CHUNK_SIZE>
class Gravitation {
	public:
	const static int pair_count = (nbod*(nbod-1))/2;

	typedef GravitationScalars<CHUNK_SIZE> shared_data [pair_count][3];

	private:
	public: // hack to make this public to test whether it's worth using shared memory for some small steps
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

	GENERIC void sum(int b,int c,double& acc, double & jerk)const
        { sum_works(b,c,acc,jerk); }

	/*  
	 *  Find the acceleration and jerk for a planet.
	 *
	 *  @b  planet number
	 *  @c  coordinate number x:0,y:1,z:2
	 *  @acc  reference to output variable for acceleration
	 *  @jerk reference to output variable for jerk
	 */
	GENERIC void sum_test(int b,int c,double& acc, double & jerk)const{
	  // Total acceleration from sun and other planets
	  double accs[2] = {0.0, 0.0};
	  // Total jerk from sun and other planets
	  double jerks[2] = {0.0, 0.0};

#pragma unroll
		for(int d = 0; d < pair_count; d++){
		  int x = first(d), y= second(d);
		  if( (x==b) || (y==b) ) 
		    {
		      double mass;
		      int dest;
		      if(x==b) 
			{
			  mass = sys[y].mass();
			  dest = (y==0) ? 0 : 1;
			}
		      else
			{
			  mass = -sys[x].mass();
			  dest = (x==0) ? 0 : 1;
			}

		      accs[dest] += shared[d][c].acc() * mass;
		      jerks[dest] += shared[d][c].jerk() * mass;
		    }
		}
		acc = accs[0] + accs[1];
		jerk = jerks[0] + jerks[1];
	}

	/*  
	 *  Find the acceleration and jerk for a planet.
	 *
	 *  @b  planet number
	 *  @c  coordinate number x:0,y:1,z:2
	 *  @acc  reference to output variable for acceleration
	 *  @jerk reference to output variable for jerk
	 */
	GENERIC void sum_works(int b,int c,double& acc, double & jerk)const{
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
				
  // TODO: Remove once allow propagators to use GravitationAccOnly
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
		 return pair_count * 3  * sizeof(GravitationScalars<CHUNK_SIZE>)/CHUNK_SIZE;
	}

	static __device__ void * system_shared_data_pointer(const int sysid_in_block) {
		extern __shared__ char shared_mem[];
		int b = sysid_in_block / CHUNK_SIZE ;
		int i = sysid_in_block % CHUNK_SIZE ;
		int idx = i * sizeof(double) 
			+ b * CHUNK_SIZE 
			* shmem_per_system();
		return &shared_mem[idx];
	}

    // WARNING: Need to test that this works (accounting for larger memory usage due to coalesced arrys)
	static __device__ void * unused_shared_data_pointer(const int system_per_block) {
		extern __shared__ char shared_mem[];
		//		int idx = system_per_block * shmem_per_system();
		int b = system_per_block / CHUNK_SIZE ;
		int i = system_per_block % CHUNK_SIZE ;
		if(i!=0) b++;
		int idx = b * CHUNK_SIZE * shmem_per_system();
		return &shared_mem[idx];
	}

				
};




/*! 
 * templatized Class to calculate acceleration and jerk in parallel
 *
 * Similar to Gravitation, but computes only terms for acceleration and not for jerk
 * To be used with integration aglorithms that don't make use of the jerk, e.g., Runge-Kutta
 */

template<int nbod, int CHUNK_SIZE = SHMEM_CHUNK_SIZE >
class GravitationAccOnly {
	public:
	const static int pair_count = (nbod*(nbod-1))/2;

	typedef GravitationAccOnlyScalars<CHUNK_SIZE> shared_data [pair_count][3];

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



};

/** 
 * Doesn't-need-to-be-templatized Class working as a function object to 
 * calculate acceleration and jerk in parallel for many-body systems.
 * 
 * shared_data member is to be initialized with a shared memory pointer.
 * Although it does not have to be a shared memory pointer.
 *
 * It operates in two steps:
 *
 * Step 0: Write positions and velocities to global memory which is cached. 
 *
 * Step 1: Loops over other bodies to: 
 * Step 1a: Calculate distances between pairs of bodies
 * The intermediate values for calculating acc/jerk are stored to the
 * shared memory.
 *
 * Step 1b: Calculate forces that affect the thread's body using data 

 * This function should be called per body per component. It uses the shared
 * data and calculates the acc/jerk for each body.
 *
 */
template<int nbod, int CHUNK_SIZE = SHMEM_CHUNK_SIZE>
class GravitationLargeN {
	public:
        const static int body_count = nbod;

        typedef GravitationScalars<CHUNK_SIZE> shared_data [body_count];

	private:
	public: // hack to make this public to test whether it's worth using shared memory for some small steps
	ensemble::SystemRef& sys;
	shared_data &shared;

	/**
	 * Helper function for calculating inner product
	 */
	GENERIC static double inner_product(const double a[3],const double b[3]){
		return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
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
	GENERIC GravitationLargeN(ensemble::SystemRef& sys,shared_data &shared):sys(sys),shared(shared){	}

	/*  
	 *  Find the acceleration and jerk for a planet.
	 *
	 *  @b  planet number
	 *  @c  coordinate number x:0,y:1,z:2
	 *  @acc  reference to output variable for acceleration
	 *  @jerk reference to output variable for jerk
	 */
         GPUAPI void sum(int b,int c,double& acc, double & jerk)const{

	   // Total acceleration 
	   acc = 0.0;
	   // Total jerk 
	   jerk = 0.0;

	   // Should we try to load these into registers? 
	   // Or is it hopeless/cache good enough anyway
	   //	   double pos_this_body = sys[b][c].pos();
	   //	   double vel_this_body = sys[b][c].vel();

	   // Loop over all other bodies, saving sun for last
	   for(int bb=nbod-1;bb>=0;--bb)
	     {
	       if( (sys[bb].mass()>0.0) && (b!=bb) && (c==0) )
		 {
		   // Relative vector from planet b to planet bb
		   double dx[3] =  { sys[bb][0].pos()- sys[b][0].pos(),sys[bb][1].pos()- sys[b][1].pos(), sys[bb][2].pos()- sys[b][2].pos() };
		   // Relative velocity between the planets
		   double dv[3] =  { sys[bb][0].vel()- sys[b][0].vel(),sys[bb][1].vel()- sys[b][1].vel(), sys[bb][2].vel()- sys[b][2].vel() };
		   // Distance between the planets
		   double r2 =  dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2];
		   
		   // broadcast jerk_mag and acc_mag via shared memory
		   shared[b].jerk() = inner_product(dx,dv) * 3.0 / r2;
		   shared[b].acc() = sys[bb].mass() * rsqrt(r2) / r2;
		 
		 }
	       __syncthreads();
	       if( (sys[bb].mass()>0.0) && (b!=bb) )
		 {
		   // Relative vector from planet bb to planet b
		   double dx =  sys[bb][c].pos()- sys[b][c].pos();
		   // Relative velocity between the planets
		   double dv = sys[bb][c].vel()- sys[b][c].vel();
		   acc += dx* shared[b].acc();
		   jerk += (dv - dx * shared[b].jerk() ) * shared[b].acc();
		 }
	     }

	 }


	/*  
	 *  Find the acceleration for a planet.
	 *
	 *  @b  planet number
	 *  @c  coordinate number x:0,y:1,z:2
	 *
	 * \todo Remove once allow propagators to use GravitationAccOnly
	 *
	 */
         GPUAPI double sum_acc(int b,int c)const{


	   // Total acceleration 
	   double acc = 0.0;
	   // Total jerk 
	   //	   jerk = 0.0;

	   // Should we try to load these into registers? 
	   // Or is it hopeless/cache good enough anyway
	   //	   double pos_this_body = sys[b][c].pos();
	   //	   double vel_this_body = sys[b][c].vel();

	   // Loop over all other bodies, saving sun for last
	   for(int bb=nbod-1;bb>=0;--bb)
	     {
	       if( (sys[bb].mass()>0.0) && (b!=bb) && (c==0) )
		 {
		   // Relative vector from planet b to planet bb
		   double dx[3] =  { sys[bb][0].pos()- sys[b][0].pos(),sys[bb][1].pos()- sys[b][1].pos(), sys[bb][2].pos()- sys[b][2].pos() };
		   // Relative velocity between the planets
		   double dv[3] =  { sys[bb][0].vel()- sys[b][0].vel(),sys[bb][1].vel()- sys[b][1].vel(), sys[bb][2].vel()- sys[b][2].vel() };
		   // Distance between the planets
		   double r2 =  dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2];
		   
		   // broadcast jerk_mag and acc_mag via shared memory
		   //		   shared[b].jerk() = inner_product(dx,dv) * 3.0 / r2;
		   shared[b].acc() = sys[bb].mass() * rsqrt(r2) / r2;
		 
		 }
	       __syncthreads();
	       if( (sys[bb].mass()>0.0) && (b!=bb) )
		 {
		   // Relative vector from planet bb to planet b
		   double dx =  sys[bb][c].pos()- sys[b][c].pos();
		   // Relative velocity between the planets
		   double dv = sys[bb][c].vel()- sys[b][c].vel();
		   acc += dx* shared[b].acc();
		   //		   jerk += (dv - dx * shared[b].jerk() ) * shared[b].acc();
		 }
	     }
	   return acc;

	 }


	/*  
	 *  Find the acceleration for a planet, 
	 *  ignoring body 0 (presumably central star)
	 *
	 *  @b  planet number
	 *  @c  coordinate number x:0,y:1,z:2
	 *
	 * \todo Remove once allow propagators to use GravitationAccOnly
	 *
	 */

         GPUAPI double sum_acc_planets(int b,int c)const{

	   // Total acceleration 
	   double acc = 0.0;
	   // Total jerk 
	   //	   jerk = 0.0;

	   // Should we try to load these into registers? 
	   // Or is it hopeless/cache good enough anyway
	   //	   double pos_this_body = sys[b][c].pos();
	   //	   double vel_this_body = sys[b][c].vel();

	   // Loop over all other bodies, saving sun for last
	   for(int bb=nbod-1;bb>=1;--bb)
	     {
	       if( (sys[bb].mass()>0.0) && (b!=bb) && (c==0) )
		 {
		   // Relative vector from planet b to planet bb
		   double dx[3] =  { sys[bb][0].pos()- sys[b][0].pos(),sys[bb][1].pos()- sys[b][1].pos(), sys[bb][2].pos()- sys[b][2].pos() };
		   // Relative velocity between the planets
		   double dv[3] =  { sys[bb][0].vel()- sys[b][0].vel(),sys[bb][1].vel()- sys[b][1].vel(), sys[bb][2].vel()- sys[b][2].vel() };
		   // Distance between the planets
		   double r2 =  dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2];
		   
		   // broadcast jerk_mag and acc_mag via shared memory
		   //		   shared[b].jerk() = inner_product(dx,dv) * 3.0 / r2;
		   shared[b].acc() = sys[bb].mass() * rsqrt(r2) / r2;
		 
		 }
	       __syncthreads();
	       if( (sys[bb].mass()>0.0) && (b!=bb) )
		 {
		   // Relative vector from planet bb to planet b
		   double dx =  sys[bb][c].pos()- sys[b][c].pos();
		   // Relative velocity between the planets
		   double dv = sys[bb][c].vel()- sys[b][c].vel();
		   acc += dx* shared[b].acc();
		   //		   jerk += (dv - dx * shared[b].jerk() ) * shared[b].acc();
		 }
	     }
	   return acc;
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
		if(b < nbod && c < 3){
			sum(b,c,acc,jerk);
		}
	}

  // TODO: Remove once allow propagators to use GravitationAccOnly
	__device__ double acc (int ij,int b,int c,double& pos,double& vel)const{
		// Write positions to shared (global) memory
		if(b < nbod && c < 3)
			sys[b][c].pos() = pos , sys[b][c].vel() = vel;
		__syncthreads();
		if(b < nbod && c < 3){
		  return sum_acc(b,c);
		}
		else { return 0; }
	}

  // TODO: Remove once allow propagators to use GravitationAccOnly
	__device__ double acc_planets (int ij,int b,int c,double& pos,double& vel)const{
		// Write positions to shared (global) memory
		if(b < nbod && c < 3)
			sys[b][c].pos() = pos , sys[b][c].vel() = vel;
		__syncthreads();
		if(b < nbod && c < 3){
		  return sum_acc_planets(b,c);
		}
		else { return 0; }
	}


	static GENERIC int shmem_per_system() {
	  const int body_count = nbod;
	  return body_count * sizeof(GravitationScalars<CHUNK_SIZE>)/CHUNK_SIZE;
	}

	static __device__ void * system_shared_data_pointer(const int sysid_in_block) {
		extern __shared__ char shared_mem[];
		int b = sysid_in_block / CHUNK_SIZE ;
		int i = sysid_in_block % CHUNK_SIZE ;
		int idx = i * sizeof(double) 
			+ b * CHUNK_SIZE 
			* shmem_per_system();
		return &shared_mem[idx];
	}

    // WARNING: Need to test that this works (accounting for larger memory usage due to coalesced arrys)
	static __device__ void * unused_shared_data_pointer(const int system_per_block) {
		extern __shared__ char shared_mem[];
		int b = system_per_block / CHUNK_SIZE ;
		int i = system_per_block % CHUNK_SIZE ;
		if(i!=0) b++;
		int idx = b * CHUNK_SIZE * shmem_per_system();
		return &shared_mem[idx];
	}

				
};



/** 
 * templatized Class working as a function object to 
 * calculate acceleration and jerk in parallel.
 *
 * \todo For this to be useful would need to make it so that integrators using this gravitation class would know that they only needed to launch 3*nbod threads per systems.  
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
template<int nbod, int CHUNK_SIZE = SHMEM_CHUNK_SIZE>
class GravitationMediumN {
	public:
	const static int pair_count = (nbod*(nbod-1))/2;

	typedef GravitationScalars<CHUNK_SIZE> shared_data [pair_count][3];

	private:
	public: // hack to make this public to test whether it's worth using shared memory for some small steps
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
	GENERIC GravitationMediumN(ensemble::SystemRef& sys,shared_data &shared):sys(sys),shared(shared){	}

	/**
	 *  Step one of the algorithm. All pairs run in parallel. This
	 *  function calculates intermediate results for a pair and
	 *  stores it into shared array.
	 *
	 *  @arg ij Integer number refering to a pair
	 */
	GENERIC void calc_pair(int ij_first )const{
	  for(int ij = ij_first ; ij < std::max(3*nbod, pair_count) ; ij += 3*nbod )
	    {
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

	GENERIC void sum(int b,int c,double& acc, double & jerk)const
	//        { sum_test(b,c,acc,jerk); }
	        { sum_works(b,c,acc,jerk); }

	/*  
	 *  Find the acceleration and jerk for a planet.
	 *
	 *  @b  planet number
	 *  @c  coordinate number x:0,y:1,z:2
	 *  @acc  reference to output variable for acceleration
	 *  @jerk reference to output variable for jerk
	 */
	GENERIC void sum_test(int b,int c,double& acc, double & jerk)const{
	  // Total acceleration from sun and other planets
	  double accs[2] = {0.0, 0.0};
	  // Total jerk from sun and other planets
	  double jerks[2] = {0.0, 0.0};

#pragma unroll
		for(int d = 0; d < pair_count; d++){
		  int x = first(d), y= second(d);
		  if( (x==b) || (y==b) ) 
		    {
		      double mass;
		      int dest;
		      if(x==b) 
			{
			  mass = sys[y].mass();
			  dest = (y==0) ? 0 : 1;
			}
		      else
			{
			  mass = -sys[x].mass();
			  dest = (x==0) ? 0 : 1;
			}

		      accs[dest] += shared[d][c].acc() * mass;
		      jerks[dest] += shared[d][c].jerk() * mass;
		    }
		}
		acc = accs[0] + accs[1];
		jerk = jerks[0] + jerks[1];
	}

	/*  
	 *  Find the acceleration and jerk for a planet.
	 *
	 *  @b  planet number
	 *  @c  coordinate number x:0,y:1,z:2
	 *  @acc  reference to output variable for acceleration
	 *  @jerk reference to output variable for jerk
	 */
	GENERIC void sum_works(int b,int c,double& acc, double & jerk)const{
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
	GPUAPI void operator() (int ij,int b,int c,double& pos,double& vel,double& acc,double& jerk) const{
		// Write positions to shared (global) memory
		if(b < nbod && c < 3)
			sys[b][c].pos() = pos , sys[b][c].vel() = vel;
		__syncthreads();
		
	        // if(ij < pair_count)
		if(ij < 3*nbod )
			calc_pair(ij);
		__syncthreads();
		if(b < nbod && c < 3){
			sum(b,c,acc,jerk);
		}
	}
				
  // TODO: Remove once allow propagators to use GravitationAccOnly
	__device__ double acc_planets (int ij,int b,int c)const{
	  //		if(ij < pair_count)
		if(ij < 3*nbod )
			calc_pair(ij);
		__syncthreads();
		if(b < nbod && c < 3){
			return sum_acc_planets(b,c);
		}else
			return 0;
	}
				
  // TODO: Remove once allow propagators to use GravitationAccOnly
	__device__ double acc (int ij,int b,int c,double& pos,double& vel)const{
		// Write positions to shared (global) memory
		if(b < nbod && c < 3)
			sys[b][c].pos() = pos , sys[b][c].vel() = vel;
		__syncthreads();
		//		if(ij < pair_count)
		if(ij < 3*nbod )
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
		 return pair_count * 3  * sizeof(GravitationScalars<CHUNK_SIZE>)/CHUNK_SIZE;
	}

	static __device__ void * system_shared_data_pointer(const int sysid_in_block) {
		extern __shared__ char shared_mem[];
		int b = sysid_in_block / CHUNK_SIZE ;
		int i = sysid_in_block % CHUNK_SIZE ;
		int idx = i * sizeof(double) 
			+ b * CHUNK_SIZE 
			* shmem_per_system();
		return &shared_mem[idx];
	}

    // WARNING: Need to test that this works (accounting for larger memory usage due to coalesced arrys)
	static __device__ void * unused_shared_data_pointer(const int system_per_block) {
		extern __shared__ char shared_mem[];
		//		int idx = system_per_block * shmem_per_system();
		int b = system_per_block / CHUNK_SIZE ;
		int i = system_per_block % CHUNK_SIZE ;
		if(i!=0) b++;
		int idx = b * CHUNK_SIZE * shmem_per_system();
		return &shared_mem[idx];
	}

				
};





}
}
}

