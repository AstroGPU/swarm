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

/*! \file log.h
 *  \brief declares hlog & dlog
 *
*/

#pragma once

#include "gpulog/gpulog.h"
#include "gpulog/lprintf.h"

#include "swarm.h"

namespace swarm
{

	// for on-GPU state logging of bodies
	// TODO: Move this to swarm.h
	// NOTE: I've written out the datatypes _explicitly_, because
	// of alignment requirements that have to be hand-tuned between
	// the device and host code. Yes, this _is_ unfortunate.
	struct ALIGN(8) body
	{
		// NOTE: put all doubles first, to avoid interstitial padding
		// and alignment nvcc vs. gcc issues
		double  x, y, z;
		double  vx, vy, vz;
		double   mass;
		long bod;

		// load body information from ensemble to body structure
		GENERIC void set(const int& i,const ensemble::Body& b)
		{
			bod = i;
			mass = b.mass();
			x = b[0].pos();
			y = b[1].pos();
			z = b[2].pos();
			vx = b[0].vel();
			vy = b[1].vel();
			vz = b[2].vel();
		}
	};

	// body_set class: hold a set of indices to bodies in a given system in
	// a given ensemble. This class should be constructed using make_body_set() functions and
	// used in conjunction with log::event to store a set of bodies in a single event
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

        /// swarm logging system.  See docs/eventlog.html
	namespace log
	{
		static const int EVT_SNAPSHOT		= 1;	// marks a snapshot of a system. see swarm::log::system() down below
		static const int EVT_EJECTION		= 2;	// marks an ejection event


		template<typename L, typename T1>
GENERIC PTR_T(SCALAR(T1)) event(L &l, const int recid, const double T, const int sys, const T1 &v1)
			{
				return l.write(recid, T, sys, v1);
			}
	
		template<typename L, typename T1, typename T2>
GENERIC PTR_T(SCALAR(T2)) event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2)
			{
				return l.write(recid, T, sys, v1, v2);
			}
	
		template<typename L, typename T1, typename T2, typename T3>
GENERIC PTR_T(SCALAR(T3)) event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2, const T3 &v3)
			{
				return l.write(recid, T, sys, v1, v2, v3);
			}
	
		template<typename L, typename T1, typename T2, typename T3, typename T4>
GENERIC PTR_T(SCALAR(T4)) event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4)
			{
				return l.write(recid, T, sys, v1, v2, v3, v4);
			}
	
		template<typename L, typename T1, typename T2, typename T3, typename T4, typename T5>
GENERIC PTR_T(SCALAR(T5)) event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5)
			{
				return l.write(recid, T, sys, v1, v2, v3, v4, v5);
			}
	
		template<typename L, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
GENERIC PTR_T(SCALAR(T6)) event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6)
			{
				return l.write(recid, T, sys, v1, v2, v3, v4, v5, v6);
			}
	
		template<typename L, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
GENERIC PTR_T(SCALAR(T7)) event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7)
			{
				return l.write(recid, T, sys, v1, v2, v3, v4, v5, v6, v7);
			}
	
		template<typename L, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
GENERIC PTR_T(SCALAR(T8)) event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7, const T8 &v8)
			{
				return l.write(recid, T, sys, v1, v2, v3, v4, v5, v6, v7, v8);
			}


		/*
			Store a snapshot of the entire system (EVT_SNAPSHOT).
		*/
		template<typename L>
		GENERIC void system(L &l, ensemble::SystemRefConst sys)
		{
			body *bodies = swarm::log::event(l, EVT_SNAPSHOT, sys.time(), sys.number() , sys.flags() , sys.nbod(), gpulog::array<body>(sys.nbod()));
			if(bodies != NULL) // buffer overflow hasn't happened
			{
				for(int bod=0; bod != sys.nbod(); bod++)
				{
					bodies[bod].set(bod,sys[bod]);
				}
			}
		}

		/*
			Store a snapshot of the entire system (EVT_SNAPSHOT).
		*/
		template<typename L>
		GENERIC void system(L &l, const ensemble &ens, const int sys)
		{
			system(l,ens[sys]);
		}

		/*
			Store a snapshot of the entire ensemble (convenience).
		*/
		template<typename L>
		GENERIC void ensemble(L &l, const swarm::ensemble &ens)
		{
			for(int sys = 0; sys < ens.nsys(); sys++)
			{
				system(l, ens[sys]);
			}
		}

	}
}

/// For use by gpu logging subsystem.  See docs/eventlog.html
namespace gpulog
{
	namespace internal
	{
	 	// body_set_cls is a proxy for an array of bodies, so make sure it reports
		// the same alignment as body[N], as well as sizeof()
		template<int N> struct alignment<swarm::body_set<N> > : public alignment<swarm::body[N]> { };	// return alignment of body[N]
		template<int N> struct    ttrait<swarm::body_set<N> > : public    ttrait<swarm::body[N]> { };	// return traits of body[N]

		// serialization of a list of bodies
		template<int N> struct     argio<swarm::body_set<N> >
		{
			GENERIC static void put(char *ptr, const swarm::body_set<N> &br, int start, int datalen)
			{
				DHOST( std::cerr << "Writing [" << br << "] start=" << start << " len=" << datalen << "\n" );
				DGPU( printf("Writing start=%d len=%d\n", start, datalen); );
				dev_assert(sizeof(swarm::body)*N == datalen);

				// write out N bodies
				swarm::body *bodies = (swarm::body *)(ptr + start);
				for(int i=0; i != N; i++)
				{
					bodies[i].set(i,br.ens[br.sys][br.bod[i]]);
				}
			}
		};
	}
}

