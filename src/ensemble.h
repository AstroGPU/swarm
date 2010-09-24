/*************************************************************************
 * Copyright (C) 2010 by Mario Juric  and the Swarm-NG Development Team  *
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

/*! \file ensemble.h
 *   \brief Declares ensemble class 
 *
*/

#include <stdexcept>
#include <string>
#include <cstring>
#include <map>
#include <cassert>
#include <cmath>
#include <vector>

#ifdef THROW_IS_ABORT
	#include <cassert>
	#include <cstring>
        #include <cstdio>
#endif

/// The main namespace for the Swarm-NG library
namespace swarm {

class integrator;
class gpu_ensemble;
class cpu_ensemble;

typedef double real_time;
typedef float  real_mass;
typedef double real_pos;
typedef double real_vel;
typedef unsigned int uint;

/**
	\brief Ensemble data storage base class

	Ensemble base class. Defines the minimum subset of data commonly needed
	by all integrators (as well as the user). Must not be used directly by
	the user; use cpu_ensemble and gpu_ensemble instead.

	Note that it has no constructors/destructors/virtuals, to allow its
	instantiation in __constant__ memory on the GPU.
*/

class ensemble
{
	public:
		enum { INACTIVE = 0x01 };

	protected:
		int m_nsys, m_nbod; //!< number of active (currently integrating) systems, total number of systems and number of bodies per system

		// m_nsys wide array
		real_time *m_T;
		real_time *m_Tend;
		uint *m_nstep;

		// m_nsys*2 wide arrays
		real_time *m_Toutput;

		// m_nsys*m_nbod*3 wide arrays
		real_pos  *m_xyz;
		real_vel  *m_vxyz;

		// m_nsys*m_nbod wide arrays
		real_mass *m_m;

		// m_nsys wide arrays
		int	*m_flags;			//!< flags about a given system. Bit 0 is the inactivity flag (m_flags & 0x01 ? inactive : acctive). Others can be user-defined.
		int	*m_systemIndices;		//!< map from original systemID to sys index in m_* arrays
		integrator *m_last_integrator;

		// scalars
		//int	*m_nactive;			//!< number of active systems, sum of !(m_flags & 0x01). Computed on GPU on exit from kernel.

	public:
		// DEPRECATED: For Young In's code ONLY!!!
		real_pos *xyz()  { return m_xyz; }
		real_vel *vxyz() { return m_vxyz; }
		const real_pos *xyz()  const { return m_xyz; }
		const real_vel *vxyz() const { return m_vxyz; }

	protected:
		void construct_base()
		{
			memset(this, 0, sizeof(*this));
		}

	public:


		//! Point structure delegating items from the ensemble
		template<class real_pos>
		struct point {
			real_pos x,y,z;
			__device__ point(){
			};
			__host__ __device__ point(const real_pos x ,const real_pos y ,const real_pos z ):x(x),y(y),z(z){};
			template<class T>
			__device__ point (const point<T>& a):x(a.x),y(a.y),z(a.z) {}
			__device__ real_pos * data() {
				return (real_pos *) this;
			}
			__device__ const real_pos * data()const {
				return (real_pos *) this;
			}
			template<class T>
			__device__ point operator -(const point<T>& a)const {
				return point(x-a.x,y-a.y,z-a.z);
			}
			template<class T>
			__device__ point operator +(const point<T>& a)const {
				return point(x+a.x,y+a.y,z+a.z);
			}
			template<class T>
			__device__ point operator *(const T& a)const {
				return point(x*a,y*a,z*a);
			}
			template<class T>
			__device__ point operator /(const T& a)const {
				return point(x/a,y/a,z/a);
			}
			template<class T>
			__device__ point& operator +=(const point<T>& a){
				x += a.x, y += a.y, z+= a.z;
				return (*this);
			}
			template<class T>
			__device__ point& operator *=(const T& a){
				x *= a, y *= a, z*= a;
				return (*this);
			}
			template<class T>
			__device__ point& operator -=(const point<T>& a){
				x -= a.x, y -= a.y, z-= a.z;
				return (*this);
			}
			__device__ real_pos length() const {
				return sqrt(x*x + y*y + z*z);
			}
			__device__ real_pos squared_length() const {
				return (x*x + y*y + z*z);
			}
			__device__ point xy_project()const{
				return point(x,y,0);
			}
			__device__ point normalize()const{
				return (*this)/length();
			}
			__device__ point outer_product(const point& a)const {
				return point( y*a.z - a.y*z, z*a.x - x*a.z, x*a.y-y*a.x);
			}
			__device__ real_pos operator *(const point& a)const{
				return x*a.x + y*a.y + z*a.z;
			}
			template<class T>
			__device__ point& operator =(const point<T>& a){
				x = a.x, y = a.y, z = a.z;
				return (*this);
			}
			//static point i,j,k,o;
		};
		/*
		struct pointref {
			real_pos &x,&y,&z;
			__device__ pointref(real_pos &xx,real_pos& yy,real_pos &zz):x(xx),y(yy),z(zz){};
			__device__ pointref(const pointref& p):x(p.x),y(p.y),z(p.z){}
			__device__ operator point<real_pos>()const{ return point<real_pos>(x,y,z);}
			__device__ pointref& operator =(const point<real_pos>& p){ x = p.x, y = p.y, z = p.z; return *this;}
			__device__ pointref& operator =(const pointref& p){ x = p.x, y = p.y, z = p.z; return *this;}
		};
		*/
		typedef point<real_pos> pointref;

		struct bodyref {
			ensemble * ens;
			int sys,bod;
			__host__ __device__ bodyref(ensemble& ens,const int& sys,const int& bod):ens(&ens),sys(sys),bod(bod){};
			__host__ __device__ real_pos& p(int c) { return ens->p(sys,bod,c) ; }
			__host__ __device__ real_pos& v(int c) { return ens->v(sys,bod,c) ; }
			__host__ __device__ pointref  pos() { return pointref(ens->x(sys,bod),ens->y(sys,bod),ens->z(sys,bod)) ; }
			__host__ __device__ pointref  vel() { return pointref(ens->vx(sys,bod),ens->vy(sys,bod),ens->vz(sys,bod)) ; }
			__host__ __device__ void set_pos(const pointref& p){ ens->x(sys,bod) = p.x, ens->y(sys,bod) = p.y, ens->z(sys,bod) = p.z; }
			__host__ __device__ void set_vel(const pointref& v){ ens->vx(sys,bod) = v.x, ens->vy(sys,bod) = v.y, ens->vz(sys,bod) = v.z; }
			__host__ __device__ real_pos mass() { return ens->mass(sys,bod); }
		};
		struct systemref {
			ensemble * ens;
			int sys;
			__host__ __device__ systemref(ensemble& ens,const int& sys):ens(&ens),sys(sys){};
			__host__ __device__ bodyref operator [](const int& i){ return bodyref(*ens,sys,i); }
			__host__ __device__ real_time& time() { return ens->time(sys); }
			__host__ __device__ real_time& time_end() { return ens->time_end(sys); }
			__host__ __device__ void increase_stepcount() { ens->nstep(sys)++; }
			__host__ __device__ int nbod() { return ens->nbod(); }
			/*
			__host__ __device__ pointref  pos(int bod) { return pointref(ens->x(sys,bod),ens->y(sys,bod),ens->z(sys,bod)) ; }
			__host__ __device__ pointref  vel(int bod) { return pointref(ens->vx(sys,bod),ens->vy(sys,bod),ens->vz(sys,bod)) ; }
			__host__ __device__ real_pos mass(int bod) { return ens->mass(sys,bod); }
			*/
			__host__ __device__ void set_time(const real_time& t) { ens->time(sys) = t; }
		};

		/*
		__host__ __device__ pointref  pos(int sys, int bod) { return pointref(x(sys,bod),y(sys,bod),z(sys,bod)) ; }
		*/
		__host__ __device__ systemref  sys(int sys) { return systemref(*this,sys); }
		__host__ __device__ systemref operator [](const int & i){ return sys(i); };




		// non-const versions

		/// return system time
		__host__ __device__ real_time&   time(int sys) { return m_T[sys]; }
		/// return system destination time
		__host__ __device__ real_time&   time_end(int sys) { return m_Tend[sys]; }
		/// return system ... time 
		__host__ __device__ real_time&   time_output(int sys, int k) { return m_Toutput[k*m_nsys + sys]; }
	
		/// return the current position compononte of the body  
		__host__ __device__ real_pos&  p(int sys, int bod, int c) { return m_xyz[m_nbod*m_nsys*c + bod*m_nsys + sys]; }
		/// return the current position x of the body  
		__host__ __device__ real_pos&  x(int sys, int bod) { return m_xyz[bod*m_nsys + sys]; }
		/// return the current position y of the body  
		__host__ __device__ real_pos&  y(int sys, int bod) { return m_xyz[m_nbod*m_nsys + bod*m_nsys + sys]; }
		/// return the current position z of the body  
		__host__ __device__ real_pos&  z(int sys, int bod) { return m_xyz[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }

		/// return the current velocity component of the body  
		__host__ __device__ real_vel& v(int sys, int bod, int c) { return m_vxyz[m_nbod*m_nsys*c + bod*m_nsys + sys]; }
		/// return the current velocity x of the body  
		__host__ __device__ real_vel& vx(int sys, int bod) { return m_vxyz[bod*m_nsys + sys]; }
		/// return the current velocity y of the body  
		__host__ __device__ real_vel& vy(int sys, int bod) { return m_vxyz[m_nbod*m_nsys + bod*m_nsys + sys]; }
		/// return the current velocity z of the body  
		__host__ __device__ real_vel& vz(int sys, int bod) { return m_vxyz[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }

		/// return the mass of the body  
		__host__ __device__ float& mass(int sys, int bod)   { return m_m[bod*m_nsys + sys]; }

		/// return the flags of the system  
		__host__ __device__ int& flags(int sys)	{ return m_flags[sys]; }

//		__host__ __device__ int& nactive() { return *m_nactive; }

		/// return the number of systems in the ensemble
		__host__ __device__ int& nsys() { return m_nsys; }
		/// return the number of bodies in a system 
		__host__ __device__ int& nbod() { return m_nbod; }

		/// return the index of the system 
		__host__ __device__ int& index_of_system(int sysId) { return m_systemIndices[sysId]; }
		/// return the step of the system 
		__host__ __device__ uint&   nstep(int sys) { return m_nstep[sys]; }



		// const versions
		
		/// return system time
		__host__ __device__ real_time time(int sys) const { return m_T[sys]; }
		/// return system destination time
		__host__ __device__ real_time time_end(int sys) const { return m_Tend[sys]; }
		/// return system ... time 
		__host__ __device__ real_time time_output(int sys, int k) const { return m_Toutput[k*m_nsys + sys]; }
	
		/// return the current position x of the body  
		__host__ __device__ real_pos  x(int sys, int bod) const { return m_xyz[bod*m_nsys + sys]; }
		/// return the current position y of the body  
		__host__ __device__ real_pos  y(int sys, int bod) const { return m_xyz[m_nbod*m_nsys + bod*m_nsys + sys]; }
		/// return the current position z of the body  
		__host__ __device__ real_pos  z(int sys, int bod) const { return m_xyz[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }

		/// return the current velocity x of the body  
		__host__ __device__ real_vel vx(int sys, int bod) const { return m_vxyz[bod*m_nsys + sys]; }
		/// return the current velocity y of the body  
		__host__ __device__ real_vel vy(int sys, int bod) const { return m_vxyz[m_nbod*m_nsys + bod*m_nsys + sys]; }
		/// return the current velocity z of the body  
		__host__ __device__ real_vel vz(int sys, int bod) const { return m_vxyz[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }

		/// return the mass of the body  
		__host__ __device__ float mass(int sys, int bod) const { return m_m[bod*m_nsys + sys]; }

		/// return the flags of the system  
		__host__ __device__ int flags(int sys)		const { return m_flags[sys]; }

//		__host__ __device__ int nactive() const { return *m_nactive; }
		
		/// return the number of systems in the ensemble
		__host__ __device__ int nsys() const { return m_nsys; }
		/// return the number of bodies in a system 
		__host__ __device__ int nbod() const { return m_nbod; }

		/// return the index of the system 
		__host__ __device__ int index_of_system(int sysId) const { return m_systemIndices[sysId]; }
		/// return the step of the system 
		__host__ __device__ uint nstep(int sys) const { return m_nstep[sys]; }


		// convenience

		/// check if the system is active 
		__host__ __device__ int is_active(int sys)		const { return !(m_flags[sys] & ensemble::INACTIVE); }
		/// check if the system is inactive 
		__host__ __device__ int is_inactive(int sys)		const { return m_flags[sys] & ensemble::INACTIVE; }
		/// set the system as active 
		__host__ __device__ void set_active(int sys)	{ m_flags[sys] = m_flags[sys] & ~ensemble::INACTIVE; }
		/// set the system as inactive 
		__host__ __device__ void set_inactive(int sys)	{  flags(sys) |= ensemble::INACTIVE; }

		/// set the body with the mass, position, and velocity 
		__host__ __device__ void set_body(int sys, int bod,  float m, real_pos x, real_pos y, real_pos z, real_vel vx, real_vel vy, real_vel vz)
		{
			int idx = bod*m_nsys + sys;

			m_m[idx]   =  m;
			m_xyz[idx] =  x; m_xyz[m_nbod*m_nsys + idx] =  y;  m_xyz[m_nbod*m_nsys*2 + idx] =  z;
			m_vxyz[idx]  = vx;  m_vxyz[m_nbod*m_nsys + idx] = vy;   m_vxyz[m_nbod*m_nsys*2 + idx] = vz;
		}

		/// return the mass, position, and velocity of the body 
		__host__ __device__ void get_body(int sys, int bod, float &m, real_pos &x, real_pos &y, real_pos &z, real_vel &vx, real_vel &vy, real_vel &vz) const
		{
			int idx = bod*m_nsys + sys;
			
			 m = m_m[idx];
			 x = m_xyz[idx];   y = m_xyz[m_nbod*m_nsys + idx];   z = m_xyz[m_nbod*m_nsys*2 + idx];
			vx = m_vxyz[idx]; vy = m_vxyz[m_nbod*m_nsys + idx]; vz = m_vxyz[m_nbod*m_nsys*2 + idx];
		}

		// utilities

		/// return barycenter of the system
                __host__ __device__ void get_barycenter(const int sys, real_pos& x, real_pos& y, real_pos& z, real_vel& vx, real_vel& vy, real_vel& vz, const int max_body_id) const 
                {
		  
                  x = 0.; y = 0.; z = 0.; vx = 0.; vy = 0.; vz = 0.;
                  double mass_sum = 0.;
                  for(int bod=0;bod<=max_body_id;++bod)
                    {
                      double m = mass(sys,bod);
                      x  += m* this->x(sys,bod);
                      y  += m* this->y(sys,bod);
                      z  += m* this->z(sys,bod);
                      vx += m* this->vx(sys,bod);
                      vy += m* this->vy(sys,bod);
                      vz += m* this->vz(sys,bod);
                      mass_sum += m;
                    }
                  x  /= mass_sum;
                  y  /= mass_sum;
                  z  /= mass_sum;
                  vx /= mass_sum;
                  vy /= mass_sum;
                  vz /= mass_sum;
                };

		/// return barycenter of the system
                __host__ __device__ void get_barycenter(const int sys, real_pos& x, real_pos& y, real_pos& z, real_vel& vx, real_vel& vy, real_vel& vz) const 
		{
		  get_barycenter(sys, x, y, z, vx, vy, vz, nbod()-1);
		}

		// Should these pay attention to active flag?

		/// set all systems time 
		__host__ __device__ void   set_time_all(const real_time tend) 
		{
		  for(int sys=0;sys<nsys();++sys)
		    time(sys) = tend;
		}
		/// set all systems destination time 
		__host__ __device__ void   set_time_end_all(const real_time tend) 
		{
		  for(int sys=0;sys<nsys();++sys)
		    time_end(sys) = tend;
		}
		/// advance all systems time 
		__host__ __device__ void   advance_time_end_all(const real_time dur) 
		{
		  for(int sys=0;sys<nsys();++sys)
		    time_end(sys) += dur;
		}
		/// ... 
		__host__ __device__ void   set_time_output_all(int k, const real_time tout) 
		{ 
		  for(int sys=0;sys<nsys();++sys)
		    time_output(sys,k) = tout;
		}

		/// return the total energy of the system  
		__host__ __device__ double calc_total_energy(int sys) const
		{
			double E = 0.;
			for (int bod1 = 0; bod1 != nbod(); bod1++)
			{
				float m1; real_pos x1[3], v1[3];
				get_body(sys, bod1, m1, x1[0], x1[1], x1[2], v1[0], v1[1], v1[2]);
				E += 0.5 * m1 * (v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);

				for (int bod2 = 0; bod2 < bod1; bod2++)
				{
					float m2; real_pos x2[3], v2[3];
					get_body(sys, bod2, m2, x2[0], x2[1], x2[2], v2[0], v2[1], v2[2]);
					real_pos dist = sqrt((x2[0] - x1[0]) * (x2[0] - x1[0]) + (x2[1] - x1[1]) * (x2[1] - x1[1]) + (x2[2] - x1[2]) * (x2[2] - x1[2]));

					E -= m1 * m2 / dist;
				}
			}
			return E;
		}

		/// calculate the total energy of all systems in the ensemble
		__host__ void calc_total_energy(double *E) const
		{
			for (int sys = 0; sys != nsys(); sys++)
			{
			  E[sys] = calc_total_energy(sys);
			}
		}

	public:
		integrator *last_integrator() { return m_last_integrator; }
		void set_last_integrator(integrator *li) { m_last_integrator = li; }

	private:
		friend class cpu_ensemble;
		friend class gpu_ensemble;
};
//end class ensemble
}
