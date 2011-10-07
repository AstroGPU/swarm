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

	__device__ double sum_acc_planets(int b,int c)const{
		double total = 0;

		/// Find the contribution from/to Sun first
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

	__device__ double sum_acc(int b,int c)const{
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

	__device__ void sum(int b,int c,double& acc, double & jerk)const{
		double total_a = 0.0;
		double from_sun_a = 0.0;
		double total_j = 0.0;
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


};

}
}
}
