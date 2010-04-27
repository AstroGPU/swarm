/*! \file swarmlog.h
 *  \brief declares hlog & dlog
 *
*/

#ifndef swarmlog_h__
#define swarmlog_h__

#include "gpulog/gpulog.h"
#include "gpulog/lprintf.h"

#include "swarm.h"

#if __CUDACC__
// The assumption is all CUDA code will be concatenated/included and compiled
// as a single source file (thus avoiding the creation of duplicate copies of 
// hlog and dlog)
gpulog::host_log hlog;
__constant__ gpulog::device_log dlog;
#endif

// declaration for g++-compiled sources
extern gpulog::host_log hlog;

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
		double	x, y, z;
		double	vx, vy, vz;

		float	m;
		int bod;

		// load body information from ensemble to body structure
		__device__ __host__ void set(const ensemble &ens, int sys, int bod_)
		{
			bod = bod_;
			m = ens.mass(sys, bod);
			x = ens.x(sys, bod);
			y = ens.y(sys, bod);
			z = ens.z(sys, bod);
			vx = ens.vx(sys, bod);
			vy = ens.vy(sys, bod);
			vz = ens.vz(sys, bod);
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

		__device__ __host__ inline body_set(const ensemble &ens_, int sys_) : ens(ens_), sys(sys_) { }
	};

	__device__ __host__ inline const body_set<1> make_body_set(const ensemble &ens, int sys, int bod0)
	{
		body_set<1> br(ens, sys);
		br.bod[0] = bod0;
		return br;
	}
	__device__ __host__ inline const body_set<2> make_body_set(const ensemble &ens, int sys, int bod0, int bod1)
	{
		body_set<2> br(ens, sys);
		br.bod[0] = bod0;
		br.bod[1] = bod1;
		return br;
	}
	__device__ __host__ inline const body_set<3> make_body_set(const ensemble &ens, int sys, int bod0, int bod1, int bod2)
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

		void init(const std::string &writer_cfg, int host_buffer_size = 10*1024*1024, int device_buffer_size = 10*1024*1024);
		void flush(int flags = memory);
		void shutdown();

		template<typename L, typename T1>
			__host__ __device__ inline PTR_T(SCALAR(T1)) event(L &l, const int recid, const double T, const int sys, const T1 &v1)
			{
				return l.write(recid, T, sys, v1);
			}
	
		template<typename L, typename T1, typename T2>
			__host__ __device__ inline PTR_T(SCALAR(T2)) event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2)
			{
				return l.write(recid, T, sys, v1, v2);
			}
	
		template<typename L, typename T1, typename T2, typename T3>
			__host__ __device__ inline PTR_T(SCALAR(T3)) event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2, const T3 &v3)
			{
				return l.write(recid, T, sys, v1, v2, v3);
			}
	
		template<typename L, typename T1, typename T2, typename T3, typename T4>
			__host__ __device__ inline PTR_T(SCALAR(T4)) event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4)
			{
				return l.write(recid, T, sys, v1, v2, v3, v4);
			}
	
		template<typename L, typename T1, typename T2, typename T3, typename T4, typename T5>
			__host__ __device__ inline PTR_T(SCALAR(T5)) event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5)
			{
				return l.write(recid, T, sys, v1, v2, v3, v4, v5);
			}
	
		template<typename L, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
			__host__ __device__ inline PTR_T(SCALAR(T6)) event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6)
			{
				return l.write(recid, T, sys, v1, v2, v3, v4, v5, v6);
			}
	
		template<typename L, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
			__host__ __device__ inline PTR_T(SCALAR(T7)) event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7)
			{
				return l.write(recid, T, sys, v1, v2, v3, v4, v5, v6, v7);
			}
	
		template<typename L, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
			__host__ __device__ inline PTR_T(SCALAR(T8)) event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7, const T8 &v8)
			{
				return l.write(recid, T, sys, v1, v2, v3, v4, v5, v6, v7, v8);
			}
		#if 0
		template<typename L, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
			__host__ __device__ inline PTR_T(SCALAR(T9)) event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7, const T8 &v8, const T9 &v9)
			{
				return l.write(recid, T, sys, v1, v2, v3, v4, v5, v6, v7, v8, v9);
			}
	
		template<typename L, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10>
			__host__ __device__ inline PTR_T(SCALAR(T10)) event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7, const T8 &v8, const T9 &v9, const T10 &v10)
			{
				return l.write(recid, T, sys, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10);
			}
		#endif

		/*
			Store a snapshot of the entire system (EVT_SNAPSHOT).
		*/
		template<typename L>
		__device__ __host__ inline void system(L &l, const ensemble &ens, const int sys, const double T)
		{
			body *bodies = swarm::log::event(l, EVT_SNAPSHOT, T, sys, ens.flags(sys), ens.nbod(), gpulog::array<body>(ens.nbod()));
			for(int bod=0; bod != ens.nbod(); bod++)
			{
				bodies[bod].set(ens, sys, bod);
			}
		}

		/*
			Store a snapshot of the entire ensemble (convenience).
		*/
		template<typename L>
		__device__ __host__ inline void ensemble(L &l, const swarm::ensemble &ens)
		{
			for(int sys = 0; sys != ens.nsys(); sys++)
			{
				system(l, ens, sys, ens.time(sys));
			}
		}

		__host__ __device__ inline bool needs_output(swarm::ensemble &ens, double T, int sys)
		{
			// simple output
			return T >= ens.time_output(sys, 0);
		}
		
		template<typename L>
		__host__ __device__ inline void output_system(L &log, swarm::ensemble &ens, double T, int sys)
		{
			// store the snapshot
			log::system(log, ens, sys, T);
		
			// set next stopping time -- find the next multiple
			// of dT closest to the current time, unless it's greater than
			// Tend, in which case set Tout = Tend.
			// If the current time is within 0.01% of a multiple of dT, 
			// set the _next_ multiple as the output time (otherwise
			// we'd have two outputs with practically equal times).
			const real_time &dT = ens.time_output(sys, 1);
			real_time &Tout = ens.time_output(sys, 0);
			Tout += ceil((T - Tout) / dT + 1e-4) * dT;
			if(Tout > ens.time_end(sys)) { Tout = ens.time_end(sys); }
		}
		
		/*!
		\brief output if needed
		
		...
		@param[out] log
		@param[in] ens
		@param[in] T
		@param[in] sys
		*/
		template<typename L>
		__host__ __device__ void output_system_if_needed(L &log, swarm::ensemble &ens, double T, int sys)
		{
			// simple output
			if(needs_output(ens, T, sys))
			{
		//		debug_hook();
		
				output_system(log, ens, T, sys);
			}
		}

		template<typename L>
		__host__ __device__ void output_systems_needing_output(L &log, swarm::ensemble &ens)
		{
			for(int sys = 0; sys != ens.nsys(); sys++)
			{
				real_time T = ens.time(sys);
				if(!needs_output(ens, T, sys)) { continue; }

				output_system(log, ens, T, sys);
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
			__host__ __device__ static inline void put(char *ptr, const swarm::body_set<N> &br, int start, int datalen)
			{
				DHOST( std::cerr << "Writing [" << br << "] start=" << start << " len=" << datalen << "\n" );
				DGPU( printf("Writing start=%d len=%d\n", start, datalen); );
				dev_assert(sizeof(swarm::body)*N == datalen);

				// write out N bodies
				swarm::body *bodies = (swarm::body *)(ptr + start);
				for(int i=0; i != N; i++)
				{
					bodies[i].set(br.ens, br.sys, br.bod[i]);
				}
			}
		};
	}
}

#endif
