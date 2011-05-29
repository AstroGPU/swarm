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

namespace swarm {
namespace hp {

inline __device__ int sysid(){
	return ((blockIdx.z * gridDim.y + blockIdx.y) * gridDim.x + blockIdx.x) * blockDim.y + threadIdx.y;
}
inline __device__ int sysid_in_block(){
	return threadIdx.y;
}
inline __device__ int thread_in_system() {
	return threadIdx.x;
}

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
template<int nbod>
class Gravitation {
	public:
	const static int pair_count = (nbod*(nbod-1))/2;
	struct shared_data {
		double acc[3][pair_count];
		double jerk[3][pair_count];
	} ;

	private:
	ensemble::systemref& sys;
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

	__device__ Gravitation(ensemble::systemref& sys,shared_data &shared):sys(sys),shared(shared){	}

	__device__ Gravitation(ensemble::systemref& sys,char * system_shmem):sys(sys)
		,shared(*( (struct shared_data*)  system_shmem ) )
	{}

	__device__ void calc_pair(int ij)const{
		int i = first( ij );
		int j = second( ij );
		if(i != j){

			double dx[3] =  { sys[j].p(0)- sys[i].p(0),sys[j].p(1)- sys[i].p(1), sys[j].p(2)- sys[i].p(2) };
			double dv[3] =  { sys[j].v(0)- sys[i].v(0),sys[j].v(1)- sys[i].v(1), sys[j].v(2)- sys[i].v(2) };
			double r2 =  dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2];
			double jerk_mag =  inner_product(dx,dv) * 3. / r2;
			double acc_mag =  rsqrt(r2) / r2;

#pragma unroll
			for(int c = 0; c < 3; c ++)
			{
				shared.acc[c][ij] = dx[c]* acc_mag;
				shared.jerk[c][ij] = (dv[c] - dx[c] * jerk_mag ) * acc_mag;
			}
		}

	}

	__device__ void calc_pair_acc(int ij)const{
		int i = first( ij );
		int j = second( ij );
		if(i != j){

			double dx[3] =  { sys[j].p(0)- sys[i].p(0),sys[j].p(1)- sys[i].p(1), sys[j].p(2)- sys[i].p(2) };
			double r2 =  dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2];
			double acc_mag =  rsqrt(r2) / r2;

#pragma unroll
			for(int c = 0; c < 3; c ++)
				shared.acc[c][ij] = dx[c]* acc_mag;
		}

	}

	__device__ double sum_values(double (&values)[3][pair_count] , int b,int c)const{
		double total = 0;
		double from_sun = 0;

		/// Find the contribution from/to Sun first
#pragma unroll
		for(int d = 0; d < pair_count; d++){
			int x = first(d), y= second(d);

			if(x == b){
				if(y == 0)
					from_sun += values[c][d]* sys[y].mass();
				else
					total += values[c][d]* sys[y].mass();
			}else if(y == b){
				if(x == 0)
					from_sun -= values[c][d]* sys[x].mass();
				else
					total -= values[c][d]* sys[x].mass();
			}
		}

		return from_sun + total;
	}



	__device__ double sum_values_no_sun(double (&values)[3][pair_count] , int b,int c)const{
		double total = 0;

		/// Find the contribution from/to Sun first
#pragma unroll
		for(int d = 0; d < pair_count; d++){
			int x = first(d), y= second(d);

			if(x == b){
				if(y != 0)
					total += values[c][d]* sys[y].mass();
			}else if(y == b){
				if(x != 0)
					total -= values[c][d]* sys[x].mass();
			}
		}

		return total;
	}


	__device__ double acc(int ij,int b,int c,double& pos,double& vel)const{
		// Write positions to shared (global) memory
		if(b < nbod)
			sys[b].p(c) = pos, sys[b].v(c) = vel;
		__syncthreads();
		if(ij < pair_count)
			calc_pair_acc(ij);
		__syncthreads();
		if(b < nbod)
			return sum_values(shared.acc,b,c);
		else
			return 0;
	}

	__device__ void operator() (int ij,int b,int c,double& pos,double& vel,double& acc,double& jerk) const{
		// Write positions to shared (global) memory
	        // TODO: Move outside this function. I don't think this really belongs here
		if(b < nbod)
			sys[b].p(c) = pos, sys[b].v(c) = vel;
		__syncthreads();
		if(ij < pair_count)
			calc_pair(ij);
		__syncthreads();
		if(b < nbod){
			acc =  sum_values(shared.acc,b,c);
			jerk = sum_values(shared.jerk,b,c);
		}
	}

	__device__ void operator() (int ij,int b,int c,double& acc,double& jerk) const{
	  // TODO: Do we really need this syncthreads/threadfence_block?
	  __syncthreads();
		if(ij < pair_count)
			calc_pair(ij);
		__syncthreads();
		if(b < nbod){
			acc =  sum_values(shared.acc,b,c);
			jerk = sum_values(shared.jerk,b,c);
		}
	}

	__device__ void calc_accel_Jerk_no_sun(int ij,int b,int c,double& acc,double& jerk) const{
	  // TODO: Do we really need this syncthreads/threadfence_block?
	  __syncthreads();
		if(ij < pair_count)
			calc_pair(ij);
		__syncthreads();
		if(b < nbod){
			acc =  sum_values_no_sun(shared.acc,b,c);
			jerk = sum_values_no_sun(shared.jerk,b,c);
		}
	}

	__device__ void calc_accel_no_sun(int ij,int b,int c,double& acc) const{
	  // TODO: Do we really need this syncthreads/threadfence_block?
	        __syncthreads();
		if(ij < pair_count)
			calc_pair_acc(ij);
		__syncthreads();
		if(b < nbod){
			acc =  sum_values_no_sun(shared.acc,b,c);
		}
	}


};

}
}
