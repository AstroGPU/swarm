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
template<class T>
class GravitationLargeN {
	public:
		const static int nbod = T::n;
        const static int body_count = nbod;
		const static int CHUNK_SIZE = SHMEM_CHUNK_SIZE;

        typedef GravitationAccJerkScalars<CHUNK_SIZE> shared_data [body_count];

	private:
	public: // hack to make this public to test whether it's worth using shared memory for some small steps
	ensemble::SystemRef& sys;
	shared_data &shared;

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
	 * \todo Remove once allow propagators to use GravitationAcc
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
	 * \todo Remove once allow propagators to use GravitationAcc
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

  // TODO: Remove once allow propagators to use GravitationAcc
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

  // TODO: Remove once allow propagators to use GravitationAcc
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
	  return body_count * sizeof(GravitationAccJerkScalars<CHUNK_SIZE>)/CHUNK_SIZE;
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

} } } // end of namespace bppt :: gpu :: swarm

