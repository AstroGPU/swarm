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

#pragma once 

#include "gpulog/gpulog.h"
#include "gpulog/lprintf.h"


namespace swarm
{
namespace hp 
{	
	// for on-GPU state logging of bodies
	// NOTE: I've written out the datatypes _explicitly_, because
	// of alignment requirements that have to be hand-tuned between
	// the device and host code. Yes, this _is_ unfortunate.
	struct ALIGN(8) body
	{
		// NOTE: put all doubles first, to avoid interstitial padding
		// and alignment nvcc vs. gcc issues
		double	m_x, m_y, m_z;
		double	m_vx, m_vy, m_vz;
		float	m_mass;
		int m_bod;

		// load body information from ensemble to body structure
		GPUAPI void set(const ensemble &ens, int sys, int bod_)
		{
			m_bod = bod_;
			m_mass = ens.mass(sys, bod_);
			m_x = ens.x(sys, bod_);
			m_y = ens.y(sys, bod_);
			m_z = ens.z(sys, bod_);
			m_vx = ens.vx(sys, bod_);
			m_vy = ens.vy(sys, bod_);
			m_vz = ens.vz(sys, bod_);
		}
		//		/*
		/// return reference to the current position x of the body  
		GPUAPI double&  x() { return m_x; };
		/// return reference to the current position y of the bm_ody  
		GPUAPI double&  y() { return m_y; };
		/// return reference to the current position z of the body  
		GPUAPI double&  z() { return m_z; };
		/// return reference to the current velocity x of the body  
		GPUAPI double& vx() { return m_vx; };
		/// return reference to the current velocity y of the body  
		GPUAPI double& vy() { return m_vy; };
		/// return reference to the current velocity z of the body  
		GPUAPI double& vz() { return m_vz; };
		/// return reference to the mass of the body  
		GPUAPI float& mass()   { return m_mass; };
		/// return reference to the id of the body  
		GPUAPI int& bod()   { return m_bod; };

	};

	// body_set class: hold a set of indices to bodies in a given system in
	// a given ensemble. This class should be constructed using make_body_set() functions and
	// used in conjunction with log::event to store a set of bodies in a single event
	template<int N>
	struct body_set
	{
		const ensemble &ens;
		int sys, bod[N];

		GPUAPI body_set(const ensemble &ens_, int sys_) : ens(ens_), sys(sys_) { }
	};

	GPUAPI const body_set<1> make_body_set(const ensemble &ens, int sys, int bod0)
	{
		body_set<1> br(ens, sys);
		br.bod[0] = bod0;
		return br;
	}
	GPUAPI const body_set<2> make_body_set(const ensemble &ens, int sys, int bod0, int bod1)
	{
		body_set<2> br(ens, sys);
		br.bod[0] = bod0;
		br.bod[1] = bod1;
		return br;
	}
	GPUAPI const body_set<3> make_body_set(const ensemble &ens, int sys, int bod0, int bod1, int bod2)
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

		enum { memory = 0x01, if_full = 0x02 };

		void init(const std::string &writer_cfg, int host_buffer_size = 50*1024*1024, int device_buffer_size = 50*1024*1024);
		void flush(int flags = memory);
		void shutdown();

		template<typename L, typename T1>
			GPUAPI PTR_T(SCALAR(T1)) event(L &l, const int recid, const double T, const int sys, const T1 &v1)
			{
				return l.write(recid, T, sys, v1);
			}
	
		template<typename L, typename T1, typename T2>
			GPUAPI PTR_T(SCALAR(T2)) event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2)
			{
				return l.write(recid, T, sys, v1, v2);
			}
	
		template<typename L, typename T1, typename T2, typename T3>
			GPUAPI PTR_T(SCALAR(T3)) event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2, const T3 &v3)
			{
				return l.write(recid, T, sys, v1, v2, v3);
			}
	
		template<typename L, typename T1, typename T2, typename T3, typename T4>
			GPUAPI PTR_T(SCALAR(T4)) event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4)
			{
				return l.write(recid, T, sys, v1, v2, v3, v4);
			}
	
		template<typename L, typename T1, typename T2, typename T3, typename T4, typename T5>
			GPUAPI PTR_T(SCALAR(T5)) event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5)
			{
				return l.write(recid, T, sys, v1, v2, v3, v4, v5);
			}
	
		template<typename L, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
			GPUAPI PTR_T(SCALAR(T6)) event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6)
			{
				return l.write(recid, T, sys, v1, v2, v3, v4, v5, v6);
			}
	
		template<typename L, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
			GPUAPI PTR_T(SCALAR(T7)) event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7)
			{
				return l.write(recid, T, sys, v1, v2, v3, v4, v5, v6, v7);
			}
	
		template<typename L, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
			GPUAPI PTR_T(SCALAR(T8)) event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7, const T8 &v8)
			{
				return l.write(recid, T, sys, v1, v2, v3, v4, v5, v6, v7, v8);
			}

		/*
			Store a snapshot of the entire system (EVT_SNAPSHOT).
		*/
		template<typename L>
		GPUAPI void system(L &l, const ensemble &ens, const int sys)
		{
			body *bodies = swarm::log::event(l, EVT_SNAPSHOT, ens.time(sys), sys, ens.flags(sys), ens.nbod(), gpulog::array<body>(ens.nbod()));
			if(bodies != NULL) // buffer overflow hasn't happened
			{
				for(int bod=0; bod != ens.nbod(); bod++)
				{
					bodies[bod].set(ens, sys, bod);
				}
			}
		}

		/*
			Store a snapshot of the entire ensemble (convenience).
		*/
		template<typename L>
		GPUAPI void ensemble(L &l, const swarm::ensemble &ens)
		{
			for(int sys = 0; sys != ens.nsys(); sys++)
			{
				system(l, ens, sys, ens.time(sys));
			}
		}

	}
}
}

/// For use by gpu logging subsystem.  See docs/eventlog.html
namespace gpulog
{
	namespace internal
	{
		using swarm::hp::body;
		using swarm::hp::body_set;
	 	// body_set_cls is a proxy for an array of bodies, so make sure it reports
		// the same alignment as body[N], as well as sizeof()
		template<int N> struct alignment<body_set<N> > : public alignment<body[N]> { };	// return alignment of body[N]
		template<int N> struct    ttrait<body_set<N> > : public    ttrait<body[N]> { };	// return traits of body[N]

		// serialization of a list of bodies
		template<int N> struct     argio<body_set<N> >
		{
GPUAPI static void put(char *ptr, const body_set<N> &br, int start, int datalen)
			{
				DHOST( std::cerr << "Writing [" << br << "] start=" << start << " len=" << datalen << "\n" );
				DGPU( printf("Writing start=%d len=%d\n", start, datalen); );
				dev_assert(sizeof(body)*N == datalen);

				// write out N bodies
				body *bodies = (body *)(ptr + start);
				for(int i=0; i != N; i++)
				{
					bodies[i].set(br.ens, br.sys, br.bod[i]);
				}
			}
		};
	}
}

