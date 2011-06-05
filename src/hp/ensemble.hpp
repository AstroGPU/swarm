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
#include "datatypes.hpp"

namespace swarm {
namespace hp {

const int ENSEMBLE_WARPSIZE = 16;


/**
 *
 * TODO: use allocator template parameter
 */
template< int _WARPSIZE>
class EnsembleBase {
	public:
	static const int WARPSIZE = _WARPSIZE;

	struct Body {
		struct Component {
			double _pos[WARPSIZE];
			double _vel[WARPSIZE];
			GPUAPI double& pos() { return _pos[0]; } 
			GPUAPI double& vel() { return _vel[0]; } 
		} component[3];

		double _mass[WARPSIZE];
		double _time[WARPSIZE];

		// Accessors 
		GPUAPI double& mass() { return _mass[0];  }
		GPUAPI double& time() { return _time[0];  }
		GPUAPI Component& operator[] (const int & i) { return component[i]; };
		GPUAPI const Component& operator[] (const int & i) const { return component[i]; };
	};



	struct SystemRef {
		const int _nbod;
		Body* _body;

		// Constructor
		GPUAPI SystemRef(const int& nbod,Body* body):_nbod(nbod),_body(body){}

		// Accessor
		GPUAPI Body& operator[](const int & i ) { return _body[i]; };
		GPUAPI double& time() { return _body[0].time(); }
		GPUAPI const int& nbod()const{ return _nbod;	}
	};

	typedef Body element_type;
	GPUAPI static size_t element_count(const int& nbod,const int& nsys){
		return (nsys + WARPSIZE) / WARPSIZE * nbod ;
	}

	protected:
	int _nbod;
	int _nsys;
	CoalescedStructArray< Body, double, WARPSIZE> _body;


	public:
	// Constructors
	GPUAPI EnsembleBase():_nbod(0),_nsys(0),_body(0,0) {};
	GPUAPI explicit EnsembleBase(const int& nbod, const int& nsys,Body* array):_nbod(nbod),_nsys(nsys),_body(array,element_count(nbod,nsys)){}
	~EnsembleBase() { 
		// TODO: We need reference counting before using this:
		// release();
	}
	// Accessors
	GPUAPI const int& nbod()const{ return _nbod;	}
	GPUAPI const int& nsys()const{ return _nsys;	}


	/** 
	 * Accessor function for first body of a system
	 * The hierarchy is like this:
	 *    Blocks ( nsys / WARPSIZE )
	 *    Bodies ( nbod )
	 *    Warps  ( nsys % WARPSIZE )
	 *
	 *    body index should come in the middle to 
	 *    provide efficient dynamic addressing.
	 *
	 */
	GPUAPI SystemRef operator[] (const int & i) { 
		const int sysinblock= i % WARPSIZE;
		const int blockid = i / WARPSIZE;
		const int idx = blockid * _nbod * WARPSIZE +  sysinblock;
		return SystemRef(_nbod,&_body[idx]) ;
	};

	GPUAPI void set_time( const int& sys, const double& time ) {
		SystemRef s = operator[] ( sys );
		s.time() = time;	
	}


	GPUAPI void set_body( const int& sys, const int& bod, const double& mass_planet
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
	
	GPUAPI void get_body(const int & sys, const int &  bod, double &  m
			, double& x, double & y, double& z
			, double&  vx, double&  vy, double&  vz) {
		SystemRef s = operator[] ( sys );

		x = s[bod][0].pos();
		y = s[bod][1].pos();
		z = s[bod][2].pos();

		vx = s[bod][0].vel();
		vy = s[bod][1].vel();
		vz = s[bod][2].vel();

		m = s[bod].mass();
	}

	GPUAPI Body* get() {
		return _body.get();
	}

	GPUAPI Body* begin() {
		return _body.begin();
	}

	GPUAPI Body* end() {
		return _body.end();
	}


	GPUAPI size_t size_in_bytes() {
		return element_count(_nbod,_nsys) * sizeof(element_type);
	}

	GPUAPI static size_t size_in_bytes(const int& nbod,const int& nsys){
		return element_count(nbod,nsys) * sizeof(element_type);
	}

	GPUAPI double calc_total_energy( int sys ) {
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

	GPUAPI void calc_total_energy(double* E){
		for (int sys = 0; sys != nsys(); sys++)
			E[sys] = calc_total_energy(sys);
	}

};

template< int W , template<class T> class _Allocator >
struct EnsembleAlloc : public EnsembleBase<W> {
	typedef EnsembleBase<W> Base;
	typedef EnsembleAlloc Self;
	typedef typename Base::Body Body;
	typedef _Allocator<Body> Allocator;

	EnsembleAlloc clone() {
		Body* b = Allocator::clone(Base::begin(),Base::end());
		return EnsembleAlloc(Base::nbod(),Base::nsys(),b);
	}

	static EnsembleAlloc create(const int& nbod, const int& nsys) {
		Body* b = Allocator::alloc( Base::element_count(nbod,nsys) );
		return EnsembleAlloc(nbod,nsys,b);
	}

	template< class Other > 
	Other cloneTo() {
		Other o = Other::create(Base::nbod(),Base::nsys());
		copyTo(o);
		return o;
	}

	template< class Other > 
	void copyTo(Other& o){
		assert(o.nsys() == Base::nsys() && o.nbod() == Base::nbod());
		alloc_copy(Allocator(),typename Other::Allocator(),Base::begin(),Base::end(),o.begin());
	}

	EnsembleAlloc(){}
	private:
	EnsembleAlloc(const int& nbod,const int& nsys, Body* b):Base(nbod,nsys,b){}
};

typedef EnsembleBase< ENSEMBLE_WARPSIZE > ensemble;

typedef EnsembleAlloc< ENSEMBLE_WARPSIZE , DefaultAllocator > defaultEnsemble;
typedef EnsembleAlloc< ENSEMBLE_WARPSIZE , DefaultAllocator > hostEnsemble;
typedef EnsembleAlloc< ENSEMBLE_WARPSIZE , DeviceAllocator > deviceEnsemble;

}
}
