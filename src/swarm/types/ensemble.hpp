/*************************************************************************
 * Copyright (C) 2011 by Saleh Dindar and the Swarm-NG Development Team  *
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

#include "allocators.hpp"
#include "coalescedstructarray.hpp"
#include "config.h"

namespace swarm {

template<class N>
GENERIC N sqr(const N& x) { return x*x; }

typedef long double_int;

/*! ensemble data structure containing nbody systems.
 *  
 *  Usage:
 *  Contains an ensemble of nsys() systems each containing nbod() bodies.
 *  use operator [] to peek into the data structure:
 *   - First [] to access a specific system
 *   - Second [] to access a specific body
 *   - Third [] to access a specific coordinate component
 *
 *  Example:
 *  \code
 *  ens[0][1].mass() = 1.0;  // Set mass of body 1 in system #0
 *  ens[2][1][0].pos() = 4;  // Set x cooridante of position of body #1 in system #2
 *  ens[3][2][1].vel() = 3;  // Set y coordinate of velocity of body #2 in system #3
 *  ens[1].time() = 3;       // Set time of system #1
 *
 *  \endcode
 *
 *  There are also compatibily accessors get_body(), set_body(), p(), v(), time(), ...
 *
 *  Two utility functions are provided for ease of use: calc_total_energy() and get_barycenter()
 *
 *  This class does not contain memory management routines and cannot
 *  be instantiated. It should be used as the \ref ensemble typedef only
 *  for parameter definition in functions and should be passed as a reference.
 *  To create an ensemble use one of \ref defaultEnsemble, \ref deviceEnsemble or \ref hostEnsemble
 *
 * 
 */
template< int _WARPSIZE>
class EnsembleBase {
	public:
	static const int WARPSIZE = _WARPSIZE;

	//! Concrete structure of Body 
	struct Body {
		struct Component {
			double _pos[WARPSIZE];
			double _vel[WARPSIZE];
			GENERIC double& pos() { return _pos[0]; } 
			GENERIC double& vel() { return _vel[0]; } 
			GENERIC const double& pos() const { return _pos[0]; } 
			GENERIC const double& vel() const { return _vel[0]; } 
		} component[3];

		double _mass[WARPSIZE];

		//! Mass of the body
		GENERIC double& mass() { return _mass[0];  }
		//! Mass of the body
		GENERIC const double& mass() const { return _mass[0];  }
		//! Index into components 0 is x, 1 is y, 2 is z
		//! Example b[0].pos() is x-coordinate of position 
		//! and b[3].vel() is z-coordinate of velocity.
		GENERIC Component& operator[] (const int & i) { return component[i]; };
		//! Index into components 0 is x, 1 is y, 2 is z
		//! Example b[0].pos() is x-coordinate of position 
		//! and b[3].vel() is z-coordinate of velocity.
		GENERIC const Component& operator[] (const int & i) const { return component[i]; };

		//! Distance of the planet to (0,0,0) 
        // TODO: Rename distance_to_origin_sq 
		GENERIC double radius() { 
			return sqr(operator[](0).pos()) 
				+ sqr(operator[](1).pos()) 
				+ sqr(operator[](2).pos());
		}
		//! Magnitude of velocity
        // TODO: Rename speed_sq
		GENERIC double speed() {
			return sqr(operator[](0).vel()) 
				+ sqr(operator[](1).vel()) 
				+ sqr(operator[](2).vel());
		}

		//! Get all position and velocities at once
		GENERIC void get(double& x,double & y, double & z
				, double & vx, double & vy, double & vz) {
			x = operator[](0).pos();
			y = operator[](1).pos();
			z = operator[](2).pos();
			vx = operator[](0).vel();
			vy = operator[](1).vel();
			vz = operator[](2).vel();
		}

	};

	//! Per system parameters: time and active.
	struct Sys {
		double _time[WARPSIZE];
		double_int _active[WARPSIZE];
		GENERIC double& time() { return _time[0];  }
		GENERIC const double& time() const { return _time[0];  }
		GENERIC double_int& active() { return _active[0];  }
		GENERIC const double_int& active()const { return _active[0];  }
	};




	//! Reference to a system within an ensemble
	//! Usage: SystemRef s = ens[i];
	struct SystemRef {
		//! Number of bodies, copied from ensemble
		const int _nbod;
		//! Serial number of the system in the ensemble
		const int _number;
		//! Pointer to the array of bodies
		Body* _body;
		//! Pointer to system parameters
		Sys* _sys;

		//! Only should be used by ensemble
		GENERIC SystemRef(const int& nbod,const int& number,Body* body,Sys* sys):_nbod(nbod),_number(number),_body(body),_sys(sys){}

		//! Access a body within the system
		GENERIC Body& operator[](const int & i ) const { return _body[i]; };
		//! Current time of the system
		GENERIC double& time() const { return _sys[0].time(); }
		//! Activity state
		//! True: still running, False: will not be integrated
		GENERIC double_int& active() const { return _sys[0].active(); }
		//! Same as active()
		GENERIC double_int& flags() const { return _sys[0].active(); }
		//! Number of bodies in this system
		GENERIC const int& nbod()const{ return _nbod;	}
		//! Serial number of the system in the ensemble
		GENERIC const int& number()const{ return _number;	}

		//! Distance between planet i and j in the system
		//! For a faster version c.f. \ref distance_squared_between(i,j)
		GENERIC double distance_between(const int& i , const int & j ) {
			return sqrt(distance_squared_between(i,j));
		}
		//! Distance squared between planet i and j in the system.
		GENERIC double distance_squared_between(const int& i , const int & j ) {
			const Body& b1 = _body[i], & b2 = _body[j];
			return sqr(b1[0].pos()-b2[0].pos())
				+ sqr(b1[1].pos()-b2[1].pos())
				+ sqr(b1[2].pos()-b2[2].pos());
		}
	};

	//! Constant encapsulation of SystemRef
	//! If the ens is constant use:
	//! SystemRefConst s = ens[i];
	struct SystemRefConst {
		SystemRef _ref;

		// Constructor
		GENERIC SystemRefConst(const SystemRef& ref):_ref(ref){}

		// Accessor
		GENERIC const Body& operator[](const int & i ) const { return _ref[i]; }
		GENERIC const double& time() { return _ref.time(); }
		GENERIC const int& number() { return _ref.number(); }
		GENERIC const double_int& active() { return _ref.active(); }
		GENERIC const double_int& flags() { return _ref.active(); }
		GENERIC const int& nbod()const{ return _ref.nbod();	}
		GENERIC double distance_squared_between(const int& i , const int & j ) { return _ref.distance_squared_between(i,j); }
		GENERIC double distance_between(const int& i , const int & j ) { return _ref.distance_between(i,j); }
	};

	//! Size of Body[] array required for an ensemble of size nbod,nsys
	GENERIC static size_t body_element_count(const int& nbod,const int& nsys){
		return (nsys + WARPSIZE) / WARPSIZE * nbod ;
	}

	//! Size of Sys[] array required for an ensemble of size nsys
	GENERIC static size_t sys_element_count(const int& nsys){
		return (nsys + WARPSIZE) / WARPSIZE ;
	}

	typedef CoalescedStructArray< Body, double, WARPSIZE> BodyArray;
	typedef CoalescedStructArray< Sys, double, WARPSIZE> SysArray;
	typedef typename BodyArray::PItem PBody;
	typedef typename SysArray::PItem PSys;

	protected:
	int _nbod;
	int _nsys;
	//! Coalesced array of Body(s)
	BodyArray _body;
	//! Coalesced array of Sys(s)
	SysArray _sys;

	public:



	public:
	//! Trivial constructor, creates an invalid ensemble
	GENERIC EnsembleBase():_nbod(0),_nsys(0),_body(0,0),_sys(0,0) {};
	//! Create an ensemble from pre-allocated body_array and sys_array arrays.
	GENERIC explicit EnsembleBase(const int& nbod, const int& nsys,PBody body_array, PSys sys_array):_nbod(nbod),_nsys(nsys),_body(body_array,body_element_count(nbod,nsys)),_sys(sys_array,sys_element_count(nsys)){}

	//! Number of bodies per system
	GENERIC const int& nbod()const{ return _nbod;	}
	//! Number of systems
	GENERIC const int& nsys()const{ return _nsys;	}
	//! Coalsed array of bodies. For INTERNAL use only
	BodyArray& bodies() { return _body; }
	//! Coalsed array of systems. For INTERNAL use only
	SysArray& systems() { return _sys; }


	/*! Index into systems array to access a system.
	 * Use: ens[i]
	 *
	 * Finds the first body of a system to set as Body* parameter
	 * The hierarchy is like this:
	 *    Blocks ( nsys / WARPSIZE )
	 *    Bodies ( nbod )
	 *    Warps  ( nsys % WARPSIZE )
	 *
	 *    body index should come in the middle to 
	 *    provide efficient dynamic addressing.
	 *
	 */
	GENERIC SystemRef operator[] (const int & i) { 
		const int sysinblock= i % WARPSIZE;
		const int blockid = i / WARPSIZE;
		const int idx = blockid * _nbod * WARPSIZE +  sysinblock;
		return SystemRef(_nbod,i,&_body[idx], &_sys[i] ) ;
	};

	GENERIC SystemRefConst operator[] (const int & i) const { 
		return SystemRefConst( const_cast<EnsembleBase*>(this)->operator[](i) ) ;
	};

	//// COMPATIBILITY ACCESSORS
	
	GENERIC double& mass(const int& sys, const int & bod){
		return operator[] ( sys )[bod].mass();
	}

	GENERIC double& p(const int& sys, const int & bod, const int& c){
		return operator[] ( sys )[bod][c].pos();
	}

	GENERIC double& v(const int& sys, const int & bod, const int& c){
		return operator[] ( sys )[bod][c].pos();
	}

	GENERIC double& x(const int& sys, const int& bod ) { return p(sys,bod,0); }
	GENERIC double& y(const int& sys, const int& bod ) { return p(sys,bod,1); }
	GENERIC double& z(const int& sys, const int& bod ) { return p(sys,bod,2); }

	GENERIC double& vx(const int& sys, const int& bod ) { return v(sys,bod,0); }
	GENERIC double& vy(const int& sys, const int& bod ) { return v(sys,bod,1); }
	GENERIC double& vz(const int& sys, const int& bod ) { return v(sys,bod,2); }

	GENERIC double& time( const int & sys ) {
		return operator[] ( sys ).time();
	}

	GENERIC long& flags(const int& sys){
		return operator[] ( sys ).active();
	}

	//// Const Versions
	
	GENERIC const double& mass(const int& sys, const int & bod)const{
		return operator[] ( sys )[bod].mass();
	}

	GENERIC const double& p(const int& sys, const int & bod, const int& c)const{
		return operator[] ( sys )[bod][c].pos();
	}

	GENERIC const double& v(const int& sys, const int & bod, const int& c)const{
		return operator[] ( sys )[bod][c].pos();
	}

	GENERIC const double& x(const int& sys, const int& bod )const { return p(sys,bod,0); }
	GENERIC const double& y(const int& sys, const int& bod )const { return p(sys,bod,1); }
	GENERIC const double& z(const int& sys, const int& bod )const { return p(sys,bod,2); }

	GENERIC const double& vx(const int& sys, const int& bod ) const { return v(sys,bod,0); }
	GENERIC const double& vy(const int& sys, const int& bod ) const { return v(sys,bod,1); }
	GENERIC const double& vz(const int& sys, const int& bod ) const { return v(sys,bod,2); }

	GENERIC const double& time( const int & sys ) const {
		return operator[] ( sys ).time();
	}

	GENERIC const long& flags(const int& sys)const {
		return operator[] ( sys ).active();
	}

	GENERIC void set_time( const int& sys, const double& time ) {
		SystemRef s = operator[] ( sys );
		s.time() = time;	
	}

	GENERIC bool is_active(const int& sys) const { 
		return operator[] ( sys ).active();
	}
	GENERIC bool is_inactive(const int& sys) const { 
		return ! operator[] ( sys ).active();
	}
	GENERIC void set_active(const int& sys) { 
		operator[] ( sys ).active() = true;
	}
	GENERIC void set_inactive(const int& sys) { 
		operator[] ( sys ).active() = false;
	}


	GENERIC void set_body( const int& sys, const int& bod, const double& mass_planet
			, const double& x, const double & y, const double& z
			, const double&  vx, const double&  vy, const double&  vz) {
		SystemRef s = operator[] ( sys );
		s[bod][0].pos()  = x ;
		s[bod][1].pos()  = y ;
		s[bod][2].pos()  = z ;

		s[bod][0].vel() = vx ;
		s[bod][1].vel() = vy ;
		s[bod][2].vel() = vz ;

		s[bod].mass() = mass_planet;

	}
	
	GENERIC void get_body(const int & sys, const int &  bod, double &  m
			, double& x, double & y, double& z
			, double&  vx, double&  vy, double&  vz) const {
		SystemRefConst s = operator[] ( sys );

		x = s[bod][0].pos();
		y = s[bod][1].pos();
		z = s[bod][2].pos();

		vx = s[bod][0].vel();
		vy = s[bod][1].vel();
		vz = s[bod][2].vel();

		m = s[bod].mass();
	}


	// Utilities
	//
	
	GENERIC void get_barycenter(const int& sys, double& x, double& y, double& z, double& vx, double& vy, double& vz, const int& max_body_id = 1000) const 
	{

		x = 0.; y = 0.; z = 0.; vx = 0.; vy = 0.; vz = 0.;
		double mass_sum = 0.;
		for(int bod=0;bod<=min(nbod()-1,max_body_id);++bod)
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


	//! Total energy (potential+kinetic) of a system
	GENERIC double calc_total_energy( int sys ) const {
		double E = 0.;
		for (int bod1 = 0; bod1 != nbod(); bod1++)
		{
			double m1; double x1[3], v1[3];
			get_body(sys, bod1, m1, x1[0], x1[1], x1[2], v1[0], v1[1], v1[2]);
			E += 0.5 * m1 * (v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);

			for (int bod2 = 0; bod2 < bod1; bod2++)
			{
				double m2; double x2[3], v2[3];
				get_body(sys, bod2, m2, x2[0], x2[1], x2[2], v2[0], v2[1], v2[2]);
				double dist = sqrt((x2[0] - x1[0]) * (x2[0] - x1[0]) + (x2[1] - x1[1]) * (x2[1] - x1[1]) + (x2[2] - x1[2]) * (x2[2] - x1[2]));

				E -= m1 * m2 / dist;
			}
		}
		return E;
	}

	GENERIC void calc_total_energy(double* E) const {
		for (int sys = 0; sys != nsys(); sys++)
			E[sys] = calc_total_energy(sys);
	}

	struct range_t {
		double average, min, max;
		range_t(const double& a,const double& m, const double& M)
			:average(a),min(m),max(M){}
	};

	//! Range of times of systems (average,min,max)
	//! Averages the time for all systems and finds min and max
	//! Useful to find the best value for destination time
	GENERIC range_t time_ranges() const {
		double time = operator[](0).time();
		double min = time;
		double max = time;
		double sum = time;
		for(int i= 1; i < nsys(); i++) {
			double time = operator[](i).time();
			if( time < min ) min = time;
			if( time > max ) max = time;
			sum += time;
		}
		return range_t(sum/nsys(),min,max);

	}

};

//! Allocator based version of ensemble containing memory management routines
template< int W , template<class T> class _Allocator >
struct EnsembleAlloc : public EnsembleBase<W> {
	typedef EnsembleBase<W> Base;
	typedef EnsembleAlloc Self;
	typedef typename Base::Body Body;
	typedef typename Base::Sys Sys;
	typedef _Allocator<Body> BodyAllocator;
	typedef _Allocator<Sys> SysAllocator;

	typedef boost::shared_ptr<Body> PBody;
	typedef boost::shared_ptr<Sys> PSys;

	//! Deep copy of this ensemble
	EnsembleAlloc clone() {
		return cloneTo<EnsembleAlloc>();
	}

	//! Create a new ensemble that can accomodate nsys systems with nbod bodies
	//! Arrays are allocated on the heap but ensemble structure is value-copied
	static EnsembleAlloc create(const int& nbod, const int& nsys) {
		PBody b ( BodyAllocator::alloc( Base::body_element_count(nbod,nsys) ), &BodyAllocator::free );
		PSys s ( SysAllocator::alloc( Base::sys_element_count(nsys) ), &SysAllocator::free );
		return EnsembleAlloc(nbod,nsys,b,s);
	}

	//! Clone to a different memory (e.g. GPU)
	template< class Other > 
	Other cloneTo() {
		Other o = Other::create(Base::nbod(),Base::nsys());
		copyTo(o);
		return o;
	}

	//! Copy to another ensemble located on a different memory
	template< class Other > 
	void copyTo(Other& o){
		assert(o.nsys() == Base::nsys() && o.nbod() == Base::nbod());
		alloc_copy(BodyAllocator(),typename Other::BodyAllocator(),Base::bodies().begin(),Base::bodies().end(),o.bodies().begin());
		alloc_copy(SysAllocator(),typename Other::SysAllocator(),Base::systems().begin(),Base::systems().end(),o.systems().begin());
	}

	EnsembleAlloc(){}
	private:
	EnsembleAlloc(const int& nbod,const int& nsys, PBody b, PSys s):
		Base(nbod,nsys,b.get(),s.get())
		,_body(b),_sys(s){}

	PSys _sys;
	PBody _body;
};

//! Base class of all ensembles
typedef EnsembleBase< ENSEMBLE_WARPSIZE > ensemble;

//! Default ensemble class for most of uses
typedef EnsembleAlloc< ENSEMBLE_WARPSIZE , DefaultAllocator > defaultEnsemble;
//! Ensemble allocated on host memory
typedef EnsembleAlloc< ENSEMBLE_WARPSIZE , DefaultAllocator > hostEnsemble;
//! Ensemble allocated on [GPU] device memory
typedef EnsembleAlloc< ENSEMBLE_WARPSIZE , DeviceAllocator > deviceEnsemble;

typedef hostEnsemble cpu_ensemble;
typedef deviceEnsemble gpu_ensemble;

}
