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


template<int W>
struct GravitationAccOnlyScalars {
	static const int WARPSIZE = W;
	typedef double scalar_t; 

	double _acc[WARPSIZE];

	// Accessors
	GENERIC double& acc() { return _acc[0];  }
};


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


/*! 
 * templatized Class to calculate acceleration and jerk in parallel
 *
 * It operates in two steps:
 *
 * Step 0: Write positions and velocities to global memory which is cached. 
 *
 * Step 1: Calculate distances between pairs using calc_pair
 * you should supply ij that is between 0 and n*(n-1)/2. It calculates the
 * inverse of distance squared between each pair of bodies.
 *
 * Step 2: Calculate forces that operate on each body.
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

        __device__ 
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

  	__device__ Gravitation(ensemble::SystemRef& sys,shared_data &shared):sys(sys),shared(shared){	}
  // TODO: Need to test this version
  __device__ Gravitation(ensemble::SystemRef& sys, const int  sysid_in_block):sys(sys),
							 shared(*( (shared_data*) system_shared_data_pointer( sysid_in_block)))
				{	  }

	__device__ void calc_pair(int ij)const{
		int i = first( ij );
		int j = second( ij );
		if(i != j){

			double dx[3] =  { sys[j][0].pos()- sys[i][0].pos(),sys[j][1].pos()- sys[i][1].pos(), sys[j][2].pos()- sys[i][2].pos() };
			double dv[3] =  { sys[j][0].vel()- sys[i][0].vel(),sys[j][1].vel()- sys[i][1].vel(), sys[j][2].vel()- sys[i][2].vel() };
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

  // TODO: Remove once allow propagators to use GravitationAccOnly
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

	__device__ void sum(int b,int c,double& acc, double & jerk)const{
		double acc_from_planets = 0.0;
		double acc_from_sun = 0.0;
		double jerk_from_planets = 0.0;
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


	__device__ void operator() (int ij,int b,int c,double& pos,double& vel,double& acc,double& jerk)const{
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


/*! 
 * templatized Class to calculate acceleration and jerk in parallel
 *
 * Similar to Gravitation, but computes only terms for acceleration and not for jerk
 * To be used with integration aglorithms that don't make use of the jerk, e.g., Runge-Kutta
 */

template<int nbod, int WARPSIZE = SHMEM_WARPSIZE>
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

	__device__ double acc_planets (int ij,int b,int c)const{
	  // WARNING: Should this function copy pos and vel to global first?
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
