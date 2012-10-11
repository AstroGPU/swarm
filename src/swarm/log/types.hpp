/*************************************************************************
 * Copyright (C) 2010 by Mario Juric and the Swarm-NG Development Team   *
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

/*! \file types.hpp 
 *    \brief Defines on-GPU state logging of bodies. 
 *
 *
 */

#pragma once
#include "../types/ensemble.hpp"

namespace swarm { /////////////////////////////////////////////////////////////////////////////

/** \brief for on-GPU state logging of bodies
 * NOTE: I've written out the datatypes _explicitly_, because
 * of alignment requirements that have to be hand-tuned between
 * the device and host code. Yes, this _is_ unfortunate.
 * NOTE: put all doubles first, to avoid interstitial padding
 * and alignment nvcc vs. gcc issues
 */
struct ALIGN(8) body
{
	double  x, y, z;
	double  vx, vy, vz;
	double   mass;
	int body_id; // EBF added

	//! load body information from ensemble to body structure
	GENERIC void set(const int& i,const ensemble::Body& b)
	{
		mass = b.mass();
		x = b[0].pos();
		y = b[1].pos();
		z = b[2].pos();
		vx = b[0].vel();
		vy = b[1].vel();
		vz = b[2].vel();
		body_id = i; 
	}
};

/** body_set class: hold a set of indices to bodies in a given system in
 * a given ensemble. This class should be constructed using 
 * make_body_set() functions and
 * used in conjunction with log::event to store a set of bodies in 
 * a single event
 */
template<int N>
struct body_set
{
	const ensemble &ens;
	int sys, bod[N];

	GENERIC body_set(const ensemble &ens_, int sys_) : ens(ens_), sys(sys_) { }
};

GENERIC const body_set<1> make_body_set(const ensemble &ens, int sys, int bod0)
{
	body_set<1> br(ens, sys);
	br.bod[0] = bod0;
	return br;
}
GENERIC const body_set<2> make_body_set(const ensemble &ens, int sys, int bod0, int bod1)
{
	body_set<2> br(ens, sys);
	br.bod[0] = bod0;
	br.bod[1] = bod1;
	return br;
}
GENERIC const body_set<3> make_body_set(const ensemble &ens, int sys, int bod0, int bod1, int bod2)
{
	body_set<3> br(ens, sys);
	br.bod[0] = bod0;
	br.bod[1] = bod1;
	br.bod[2] = bod2;
	return br;
}

} ////////////////////////////////////////////////////////////////////////////////////////////////



namespace gpulog { 	namespace internal {  ////////////////////////////////////////////////////

//! body_set_cls is a proxy for an array of bodies, so make sure it reports
//! the same alignment as body[N], as well as sizeof()
//! return alignment of body[N]
template<int N> struct alignment<swarm::body_set<N> > : public alignment<swarm::body[N]> { };
//! return traits of body[N]
template<int N> struct    ttrait<swarm::body_set<N> > : public    ttrait<swarm::body[N]> { };	

/** \brief Template partial specialization of argio class from gpulog for body_set
 *
 *  argio struct is specialized for correct handling of swarm::body_set<N> structs inside
 *  gpulog subsystem.
 *
 */
template<int N> struct     argio<swarm::body_set<N> >
{
	GENERIC static void put(char *ptr, const swarm::body_set<N> &br, int start, int datalen)
	{
		DHOST( std::cerr << "Writing [" << br << "] start=" << start << " len=" << datalen << "\n" );
		DGPU( printf("Writing start=%d len=%d\n", start, datalen); );
		dev_assert(sizeof(swarm::body)*N == datalen);

		// write out N bodies
		swarm::body *bodies = (swarm::body *)(ptr + start);
		for(int i=0; i < N; i++)
		{
			bodies[i].set(i,br.ens[br.sys][br.bod[i]]);
		}
	}
};

} } ///////////////////////////////////////////////////////////////////////////////////////////

