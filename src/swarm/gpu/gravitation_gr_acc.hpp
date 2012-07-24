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

#include "gravitation_common.hpp"

namespace swarm { namespace gpu { namespace bppt {

/*! 
 * templatized Class to calculate acceleration and jerk in parallel
 *
 * Similar to Gravitation, but computes only terms for acceleration and not for jerk
 * To be used with integration aglorithms that don't make use of the jerk, e.g., Runge-Kutta
 */
template<class T>
class GravitationAcc_GR {
	public:
	const static int nbod = T::n;
	const static int pair_count = (nbod*(nbod-1))/2;
	const static int CHUNK_SIZE = SHMEM_CHUNK_SIZE;

	typedef GravitationAccScalars<CHUNK_SIZE> shared_data [pair_count][3];

	private:
	ensemble::SystemRef& sys;
	shared_data &shared;
        double c2;

	public:

        // Calculate accelerations, including approximation for weak GR
        // currently hardwired for G=M_sol=AU=1, year=2pi
        // \todo read c^2 from parameter file?  read from system attribute?
        //       set system attribute from parameter file?
	__device__ GravitationAcc_GR(ensemble::SystemRef& sys,shared_data &shared):sys(sys),shared(shared)
        {	
	  c2 = 101302340.; 
        }


	__device__ void calc_pair(int ij)const{
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

          __device__ double one_over_r(const int b1, const int b2) const
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

        // Warning untested... 
        // Based on addition to Mercury from Matthew Payne
        //  Look at Benitez & Gallardo and references there-in
        // http://www.springerlink.com/content/d8r26w6610389256/fulltext.pdf
        //                      -GM
        // Additional accel =  -----   {(4GM / r - v^2) r + 4(v.r)v}
        //                     r^3c^2
        // Is it right to assing all acceleration to planet and none to star?
	__device__ double sum_acc_gr(int b,int c)const
                {
		double acc_sum = 0;
		double one_over_r;
		// avoid extra square root, by recalculating 1/r (planet-sun) from shared mem
		for(int d = 0; d < pair_count; d++)
		  {
		    int x = first<nbod>(d), y= second<nbod>(d);
		    
		    if( ((x==b)&&(y==0)) ||  ((x==0)&&(y==b)) )
		      {
			one_over_r  = shared[d][0].acc()*shared[d][0].acc();
			one_over_r += shared[d][1].acc()*shared[d][1].acc();
			one_over_r += shared[d][2].acc()*shared[d][2].acc();
		      }
		  }
		double v2 = 0., v_dot_r = 0.;
		double dx = (sys[b][0].pos()-sys[0][0].pos());
		double dv = (sys[b][0].vel()-sys[0][0].vel());
		v2 += dx*dx;
		v_dot_r += dx*dv;
		dx = (sys[b][1].pos()-sys[0][1].pos());
		dv = (sys[b][1].vel()-sys[0][1].vel());
		v2 += dx*dx;
		v_dot_r += dx*dv;
		dx = (sys[b][2].pos()-sys[0][2].pos());
		dv = (sys[b][2].vel()-sys[0][2].vel());
		v2 += dx*dx;
		v_dot_r += dx*dv;
		// assumes G=1.
		double GM = sys[0].mass();
		double f1 = (4*GM*one_over_r-v2);
		double f2 = 4*v_dot_r;
		acc_sum = -f1*(sys[b][c].pos()-sys[0][c].pos())
		          -f2*(sys[b][c].vel()-sys[0][c].vel());
		double f0 = GM*one_over_r*one_over_r*one_over_r/c2;
		acc_sum *= f0;
		return acc_sum;
		}


                // Warning: untested and known to be incomplete...
                // Calculates acceleration on planet in barycentric frame due to J2
                // Would need to figure out how to read j2... parameter file? attribute? parameter files specifying which attribute?
                // Would need to figure out how to deal with back reaction on star before we use this
                 __device__ double sum_acc_j2(int b,int c)const
                {
		double acc_sum = 0;
		double one_over_r;
		double u2;
		// avoid extra square root, by recalculating 1/r (planet-sun) from shared mem
		for(int d = 0; d < pair_count; d++)
		  {
		    int x = first<nbod>(d), y= second<nbod>(d);
		    
		    if( ((x==b)&&(y==0)) ||  ((x==0)&&(y==b)) )
		      {
			one_over_r  = shared[d][0].acc()*shared[d][0].acc();
			one_over_r += shared[d][1].acc()*shared[d][1].acc();
			one_over_r += shared[d][2].acc()*shared[d][2].acc();
			u2 = shared[d][2].acc()*shared[d][2].acc() / one_over_r;
		      }
		  }
		double dx = (sys[b][c].pos()-sys[0][c].pos());
		double mstar = sys[0].mass();
		double mpl = sys[b].mass();
		double j2 = sys[b].attribute(1); 
		double jr2 = j2*one_over_r*one_over_r;
		double tmp2 = jr2*(7.5*u2-1.5);
		double tmp3 = jr2*3.;
#if 0		
		double jr4 = j4*one_over_r*one_over_r*one_over_r*one_over_r;
		double u4 = u2*u2;
		tmp2 += jr4*(39.375*u4 - 26.25*u2 + 1.875);
		tmp3 += jr4*(17.5*u2-7.5);
#endif
#if 0		
		double jr6 = j6*one_over_r*one_over_r*one_over_r*one_over_r*one_over_r*one_over_r;
		double u6 = u4*u2;
		tmp2 += jr6*(187.6875*u6 -216.5625*u4 +59.0625*u2 -2.1875);
		tmp3 += jr6*(86.625*u4 - 78.75*u2 + 13.125);
#endif
		if(c==2) tmp2 -= tmp3;
		double tmp1 = mstar*one_over_r*one_over_r*one_over_r;
		acc_sum += dx*tmp1*tmp2;
		// \todo Need to figure out how to get this to affect sun.
		double acc_sun_c = -mpl*acc_sum/mstar; 
		return acc_sum;
		}



	__device__ double sum_acc_planets(int b,int c)const{
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

	__device__ double sum_acc(int b,int c) const{
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

		double acc_gr = sum_acc_gr(b,c);
		return acc_from_sun + acc_from_planets + acc_gr;
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

	static GENERIC int thread_per_system(){
	  return ( 3*nbod > (nbod-1)*nbod/2 ) ? nbod*3 : (nbod-1)*nbod/2;
	}

	static GENERIC int shmem_per_system() {
		 return sizeof(shared_data)/CHUNK_SIZE;
	}


};

} } } // end of namespace bppt :: gpu :: swarm


