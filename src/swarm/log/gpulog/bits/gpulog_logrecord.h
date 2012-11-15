/***************************************************************************
 *   Copyright (C) 2010 by Mario Juric                                     *
 *   mjuric@cfa.harvard.EDU                                                *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/*! \file gpulog_logrecord.h 
 *    \brief Defines methods for processing log records. 
 *
 *
 */

#ifndef bits_gpulog_logrecord_h__
#define bits_gpulog_logrecord_h__

namespace gpulog
{
namespace internal
{

	/*
		Log unserialization infrastructure	
	*/

	//! log unserialization - reads elements from a log record 
	struct logrecord
	{
		const char *ptr;
		int at;

		#if ARGINFO
		int atarg;
		__host__ __device__ inline const arginfo &get_arginfo(int arg) const { return *((arginfo *)(ptr + hdr().infos) + arg); }
		#endif

	        //! Constructor for logrecord
		__host__ __device__ inline logrecord() : ptr(NULL), at(0) { IFARGINFO(atarg = 1); }
	        //! Constructor for logrecord
		__host__ __device__ inline logrecord(const char *ptr_) : ptr(ptr_), at(hdr().dataoffset()) { IFARGINFO(atarg = 1); }

	        //! Retrurn the header pointer
		__host__ __device__ inline const header &hdr() const { return *(header *)ptr; }
	        //! 
		__host__ __device__ inline operator bool() const { return at < hdr().len; }

	        //! 
		__host__ __device__ inline int len() const { return hdr().len; }
	        //! 
		__host__ __device__ inline int msgid() const { return hdr().msgid; }
	};

	/* workaround for CUDA 2.2 template parsing bug */
	typedef logrecord& logrecord_ref;

	template<typename T>
	__device__ void assert_arg_compatible(logrecord &in, T &v)
	{
		#if ARGINFO && !__CUDACC__
		const arginfo &ai = in.get_arginfo(in.atarg);
		typedef ttrait<T> TT;

		DHOST( std::cerr << "     Log  : " << ai << "\n"; )
		DHOST( std::cerr << "     Host : arg=" << in.atarg << " align=" << TT::align << " size=" << TT::size << " dim=" << TT::dim << " isptr=" << TT::isptr << "\n"; )

		#define AITEST(a, TT, b) \
			if(a != TT::b) \
			{ \
				std::cerr << "Assertion failed GPUTypeTraits::" #a " != HostTypeTraits::" #b << " (" << a << " != " << TT::b << ")\n"; \
				std::cerr << "     Log  : " << ai << "\n"; \
				std::cerr << "     Host : arg=" << in.atarg << " align=" << TT::align << " size=" << TT::size << " dim=" << TT::dim << " isptr=" << TT::isptr << "\n"; \
				assert(a == TT::b); \
			}

		// Check for most common issues
		// TODO: We could make this more sophisticated to check for correct extraction
		// of presized arrays, as well as the last array<>.
		AITEST( ai.align, TT, align );	// check for alignment
		AITEST(  ai.size, TT, size );	// check for scalar type size match
		#undef AITEST

		if(in.atarg+1 < in.hdr().nargs) { in.atarg++; }
		#endif
	}

	//! Unserialize a scalar
	template<typename T>
	__device__ inline logrecord_ref operator >>(logrecord &in, T &v)
	{
		if(!in) { return in; }

		// compute alignment
		int offs = ASTART(in.at, ALIGNOF(T));
		in.at = offs + sizeof(T);
		assert_arg_compatible(in, v);
		dev_assert(in.at <= in.len());	/* assert we didn't run over the size of the buffer */

		v = *(T*)(in.ptr + offs);

		DHOST( std::cerr << "reading scalar POD: offs = " << offs << " val = " << v << "\n"; )

		return in;
	}

	//! specialization for sized arrays 
	template<typename T, int N>
	__device__ inline logrecord_ref operator >>(logrecord &in, const T (&v)[N])
	{
		if(!in) { return in; }

		// compute alignment
		int offs = ASTART(in.at, ALIGNOF(T));
		in.at = offs + sizeof(T)*N;
		assert_arg_compatible(in, v);
		dev_assert(in.at <= in.len());	/* assert we didn't run over the size of the buffer */

		for(int i = 0; i != N; i++)
		{
			v[i] = ((T*)(in.ptr + offs))[i];
			DHOST( std::cerr << "reading array: offs = " << offs << " N = " << N << " val[" << i << "] = " << v[i] << "\n"; )
		}

		return in;
	}

	//! specialization for char ptr (null terminated string) 
	__device__ inline logrecord_ref operator>>(logrecord &in, char *v)
	{
		if(!in) { return in; }

		// compute alignment
		DHOST( char *v0 = v; int offs = in.at; )
		assert_arg_compatible(in, v);
		while(in.ptr[in.at])
		{
			*v = in.ptr[in.at++];
			v++;
		}
		*v = 0; in.at++;

		dev_assert(in.at <= in.len());	/* assert we didn't run over the size of the buffer */

		DHOST( std::cerr << "reading character string: offs = " << offs << " val = " << v0 << " endoff = " << in.at << "\n"; )
		return in;
	}

	//! specialization for char array (we assume it'll be a null-terminated string) 
	template<int N>
	__device__ inline logrecord_ref operator>>(logrecord &in, char v[N])
	{
		return in >> (char *)v;
	}

	//! specialization for return of pointers to non-const arrays 
	template<typename T>
	__device__ inline logrecord_ref operator >>(logrecord &in, T *&v)
	{
		// intentionally force a compile time error - you MUST only be using
		// the const T* version (below), as the internal buffer to which a
		// pointer is returned is immutable
		STATIC_ASSERTION_FAILED__Your_pointer_must_be_to_const_T_and_not_just_T______(v);
		return in;
	}

	//! specialization for return of a pointer array
	template<typename T>
	__device__ inline logrecord_ref operator >>(logrecord &in, const T *&v)
	{
		if(!in) { return in; }

		// compute alignment
		int offs = ASTART(in.at, ALIGNOF(T));
		assert_arg_compatible(in, v);
		v = (const T*)(in.ptr + offs);
		in.at = in.hdr().len;		// extraction of an array of unspecified size is always the last operation

		DHOST( std::cerr << "reading unbound array : offs = " << offs << " v = " << (void*)v << " v[0] = " << *v << "\n"; )
		DHOST( std::cerr << (char*)v - in.ptr << "\n"; )

		return in;
	}

}
}

#endif // bits_gpulog_logrecord_h__
