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

/*! \file static_accjerk.hpp
 *  \brief Templatized implementation of acc jerk.
 *
*/

/// The main namespace for the Swarm-NG library
//
#include "meta.hpp"
namespace swarm {

		template<int nbod>
		class Compute_Acc {
			ensemble::systemref& sys;
			double (&pos)[3][nbod];
			double (&acc)[3][nbod];
			public:
			__device__ Compute_Acc(ensemble::systemref sys,double (&pos)[3][nbod],double (&acc)[3][nbod])
				:sys(sys),pos(pos),acc(acc) {
				#pragma unroll
					for(int i =0; i < nbod; i++)
						acc[0][i] = acc[1][i] = acc[2][i] = 0;
			}

			__device__ static void acc_updater_component(int i,int c
					,double dx[3],double scalar,double (&acc)[3][nbod]){
				acc[c][i] += dx[c]* scalar;
			}


			__device__ void operator() (int x) const {
				const int ij = (nbod*nbod-1) - x;
				const int i = ij/nbod, j = ij%nbod;
				compute_pair(i,j);
			}
			__device__ void compute_pair(int i,int j) const{
				if(i < j){

					double dx[3] =  { pos[0][j]-pos[0][i],pos[1][j]-pos[1][i],pos[2][j]-pos[2][i]};

					// computing scalar part of the acceleration
					double r2 =  dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2] ;
					double rinv = rsqrt(r2)  / r2;

					// vectorized part
					const double scalar_i = +rinv*sys[j].mass();
					acc_updater_component(i,0,dx,scalar_i,acc);
					acc_updater_component(i,1,dx,scalar_i,acc);
					acc_updater_component(i,2,dx,scalar_i,acc);

					const double scalar_j = -rinv*sys[i].mass();
					acc_updater_component(j,0,dx,scalar_j,acc);
					acc_updater_component(j,1,dx,scalar_j,acc);
					acc_updater_component(j,2,dx,scalar_j,acc);
				}
			}

			__device__ void compute(){
				//#pragma unroll
				//for(int x = 0; x < nbod*nbod ;x++) (*this)(x);
				Unroller<0,nbod*nbod>::step(*this);
			}
		};

		template<int nbod>
		class Compute_AccJerk {
			ensemble::systemref& sys;
			double (&pos)[3][nbod];
			double (&vel)[3][nbod];
			double (&acc)[3][nbod];
			double (&jerk)[3][nbod];
			public:
			__device__ Compute_AccJerk(ensemble::systemref sys
					,double (&pos)[3][nbod],double (&vel)[3][nbod]
					,double (&acc)[3][nbod],double (&jerk)[3][nbod])
				:sys(sys),pos(pos),vel(vel),acc(acc),jerk(jerk) {
				#pragma unroll
					for(int i =0; i < nbod; i++)
						acc[0][i] = acc[1][i] = acc[2][i] = 0,
						jerk[0][i] = jerk[1][i] = jerk[2][i] = 0;
			}

			__device__ static void accjerk_updater_component(int i,int c
					,double dx[3],double dv[3],double scalar,double rv
					,double (&acc)[3][nbod],double (&jerk)[3][nbod]){
				acc[c][i] += dx[c]* scalar;
				jerk[c][i] += (dv[c] - dx[c] * rv) * scalar;

			}

			__device__ void operator() (int x) const {
				const int ij = (nbod*nbod-1) - x;
				const int i = ij/nbod, j = ij%nbod;
				compute_pair(i,j);
			}

			
			inline __device__ static double inner_product(const double a[3],const double b[3]){
				return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
			}

			__device__ void compute_pair(int i,int j) const {
				if(i < j){
					double dx[3] =  { pos[0][j]-pos[0][i],pos[1][j]-pos[1][i],pos[2][j]-pos[2][i]};
					double dv[3] =  { vel[0][j]-vel[0][i],vel[1][j]-vel[1][i],vel[2][j]-vel[2][i]};

					// computing scalar part of the acceleration
					double r2 =  dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2] ;
					double rv =  inner_product(dx,dv) * 3. / r2;
					double rinv = rsqrt(r2)  / r2;

					// vectorized part
					const double scalar_i = +rinv*sys[j].mass();
					accjerk_updater_component(i,0,dx,dv,scalar_i,rv,acc,jerk);
					accjerk_updater_component(i,1,dx,dv,scalar_i,rv,acc,jerk);
					accjerk_updater_component(i,2,dx,dv,scalar_i,rv,acc,jerk);

					const double scalar_j = -rinv*sys[i].mass();
					accjerk_updater_component(j,0,dx,dv,scalar_j,rv,acc,jerk);
					accjerk_updater_component(j,1,dx,dv,scalar_j,rv,acc,jerk);
					accjerk_updater_component(j,2,dx,dv,scalar_j,rv,acc,jerk);
				}
			}

			__device__ void compute(){
				//#pragma unroll
				//for(int x = 0; x < nbod*nbod ;x++) (*this)(x);
				Unroller<0,nbod*nbod>::step(*this);
			}
		};


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
		struct shared_data {
			double acc[3][(nbod*(nbod-1))/2];
			double jerk[3][(nbod*(nbod-1))/2];
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

				// computing scalar part of the acceleration
				double r2 =  dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2];
				double rinv = rsqrt(r2)  / r2;

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

		__device__ double sum_jerk(int b,int c)const{
			double jerk = 0;
			double jerk_from_sun = 0;

			/// Find the contribution from/to Sun first
#pragma unroll
			for(int d = 0; d < (nbod*(nbod-1))/2; d++){
				int x = first(d), y= second(d);

				if(x == b){
					if(y == 0)
						jerk_from_sun += shared.jerk[c][d]* sys[y].mass();
					else
						jerk += shared.jerk[c][d]* sys[y].mass();
				}else if(y == b){
					if(x == 0)
						jerk_from_sun -= shared.jerk[c][d]* sys[x].mass();
					else
						jerk -= shared.jerk[c][d]* sys[x].mass();
				}
			}

			return jerk_from_sun + jerk;
		}

		__device__ double sum_acc(int b,int c)const{
			double	acc = 0;
			double acc_from_sun = 0;
#pragma unroll
			for(int d = 0; d < (nbod*(nbod-1))/2; d++){
				int x = first(d), y= second(d);

				if(x == b){
					if( y == 0)
						acc_from_sun  += shared.acc[c][d] * sys[y].mass();
					else
						acc  += shared.acc[c][d] * sys[y].mass();
				}else if(y == b){
					if( x == 0)
						acc_from_sun  -= shared.acc[c][d] * sys[x].mass();
					else
						acc  -= shared.acc[c][d] * sys[x].mass();
				}
			}

			return acc + acc_from_sun;
		}

		__device__ double acc(int ij,int b,int c,double& pos,double& vel)const{
			// Write positions to shared (global) memory
			if(b < nbod)
				sys[b].p(c) = pos, sys[b].v(c) = vel;
			__syncthreads();
			if(ij < (nbod*(nbod-1))/2)
				calc_pair(ij);
			__syncthreads();
			if(b < nbod)
				return sum_acc(b,c);
			else
				return 0;
		}

		__device__ void operator() (int ij,int b,int c,double& pos,double& vel,double& acc,double& jerk)const{
			// Write positions to shared (global) memory
			if(b < nbod)
				sys[b].p(c) = pos, sys[b].v(c) = vel;
			__syncthreads();
			if(ij < (nbod*(nbod-1))/2)
				calc_pair(ij);
			__syncthreads();
			if(b < nbod){
				acc = sum_acc(b,c);
				jerk = sum_jerk(b,c);
			}
		}


	};

}
