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

/*! \file gravitation_acc.hpp
 *   \brief Defines and implements class \ref swarm::gpu::bppt::GravitationAcc 
 *          that implements member functions to calculate acceleration part of the gravitation. 
 *          
 */


#pragma once

#include "gravitation_common.hpp"

namespace swarm { namespace gpu { namespace bppt {

/*! 
 * templatized Class to calculate acceleration and jerk in parallel
 *
 * Similar to Gravitation, but computes only terms for acceleration and not for jerk
 * To be used with integration aglorithms that don't make use of the jerk, e.g., Runge-Kutta
 */
template<class T>
class GravitationAcc {
	public:
	const static int nbod = T::n;
	const static int pair_count = (nbod*(nbod-1))/2;
	const static int CHUNK_SIZE = SHMEM_CHUNK_SIZE;

	typedef GravitationAccScalars<CHUNK_SIZE> shared_data [pair_count][3];

	private:
	ensemble::SystemRef& sys;
	shared_data &shared;

	public:

	GPUAPI GravitationAcc(ensemble::SystemRef& sys,shared_data &shared):sys(sys),shared(shared){	}

	private:

	GPUAPI void calc_pair(int ij)const{
		int i = first<nbod>( ij );
		int j = second<nbod>( ij );
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

	public:
          GPUAPI double one_over_r(const int b1, const int b2) const
          {
	  double sum = 0.;

		for(int d = 0; d < pair_count; d++)
		  {
		    int x = first<nbod>(d), y= second<nbod>(d);
		    if( ((x==b1) && (y==b2)) || (x==b2) && (y==b1))
		      {
			sum += shared[d][0].acc()*shared[d][0].acc();
			sum += shared[d][1].acc()*shared[d][1].acc();
			sum += shared[d][2].acc()*shared[d][2].acc();
		      }
		  }
		return sum;
	  }

	private:

	GPUAPI double sum_acc_planets(int b,int c)const{
		double acc_sum = 0;

		/// Find the contribution from/to Sun first
#pragma unroll
		for(int d = 0; d < pair_count; d++){
			int x = first<nbod>(d), y= second<nbod>(d);

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

	GPUAPI double sum_acc(int b,int c) const{
		double acc_from_planets = 0;
		double acc_from_sun = 0;

		/// Find the contribution from/to Sun first
#pragma unroll
		for(int d = 0; d < pair_count; d++){
			int x = first<nbod>(d), y= second<nbod>(d);

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

	public:

	GPUAPI void operator() (int ij,int b,int c,double& pos,double& vel,double& acc)const{
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

	/**
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

	/**
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

	static GENERIC int thread_per_system(){
	  return ( 3*nbod > (nbod-1)*nbod/2 ) ? nbod*3 : (nbod-1)*nbod/2;
	}

	static GENERIC int shmem_per_system() {
		 return sizeof(shared_data)/CHUNK_SIZE;
	}


};

} } } // end of namespace bppt :: gpu :: swarm


