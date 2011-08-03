/***************************************************************************
 *   Copyright (C) 2010 by Mario Juric   *
 *   mjuric@cfa.harvard.EDU       *
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

#ifndef bits_gpulog_ttraits_h__
#define bits_gpulog_ttraits_h__

#if GPULOG_DEBUG && !__CUDACC__
#include<iostream>
#endif

#ifdef _WIN32
#define _USE_MATH_DEFINES
#include <math.h>
#define __alignof__ __alignof
#endif 

namespace gpulog
{

	namespace internal
	{

	/* Type alignment querying */
	template<typename T> struct alignment         { static const int value = __alignof__(T); };

	/* Alignment overrides to make the host packing compatible with device packing */
	/* See table B-1 in CUDA 2.2 Programming Guide for reference */
	template<>           struct alignment<float2> { static const int value = 8; };
	template<>           struct alignment<float4> { static const int value = 16; };
	template<>           struct alignment<char2>  { static const int value = 2; };
	template<>           struct alignment<char4>  { static const int value = 4; };
	template<>           struct alignment<uchar2>  { static const int value = 2; };
	template<>           struct alignment<uchar4>  { static const int value = 4; };
	template<>           struct alignment<short2>  { static const int value = 4; };
	template<>           struct alignment<short4>  { static const int value = 8; };
	template<>           struct alignment<ushort2>  { static const int value = 4; };
	template<>           struct alignment<ushort4>  { static const int value = 8; };
	template<>           struct alignment<int2> { static const int value = 8; };
	template<>           struct alignment<int4> { static const int value = 16; };
	template<>           struct alignment<uint2> { static const int value = 8; };
	template<>           struct alignment<uint4> { static const int value = 16; };
	template<>           struct alignment<long2> { static const int value = 2*alignment<long1>::value; };	// this should take care of 32 vs. 64 bit differences, assuming -malign-double is active
	template<>           struct alignment<long4> { static const int value = 4*alignment<long1>::value; };
	template<>           struct alignment<ulong2> { static const int value = 2*alignment<ulong1>::value; };
	template<>           struct alignment<ulong4> { static const int value = 4*alignment<ulong1>::value; };
	template<>           struct alignment<double2> { static const int value = 16; };

	//
	// Type traits
	//
	template<typename T>
	struct ttrait
	{
		typedef T scalarT;

		static const size_t size = sizeof(T);
		static const size_t align = alignment<T>::value;

		static const bool isptr = false;

		static const bool isarr = false;
		
		static const bool isunspec = false;

		static const size_t dim = 1;
	};

	template<typename T>
	struct ttrait<T*> : public ttrait<T>
	{
		static const bool isptr = true;
	};

	template<typename T, int N>
	struct ttrait<T[N]> : public ttrait<T>
	{
		static const size_t dim = N;
	};

	template<typename T>
	struct ttrait<array<T> > : public ttrait<T>
	{
		static const bool isarr = true;
	};

	template<>
	struct ttrait<Tunspec>
	{
		typedef Tunspec scalarT;
		static const bool isunspec = true;

		static const size_t size = 0;
		static const size_t align = 1;
		static const bool isptr = false;
		static const bool isarr = false;
		static const size_t dim = 0;
	};

	/* Alignment, size and type trait convenience macros */
	#define ALIGNOF(T)	(ttrait<T>::align)				/* Alignment of type T */
	#define ESIZE(T)	(ttrait<T>::size)				/* Scalar element size of type T. E.g., if T = short[5], ESIZE(T) = sizeof(short) */
	#define DIMEN(T)	(ttrait<T>::dim)				/* Array size of type T. E.g., if T = short[5], DIMEN(T) = 5 */
	#define SIZEOF(T)	(ESIZE(T)*DIMEN(T))				/* Byte size of type T. For array<X>, this is the size of X */

	#define ISARRAY(T)	(ttrait<T>::isarr)				/* Is type T an array<> class */
	#define SCALAR(T)	typename gpulog::internal::ttrait<T>::scalarT	/* The scalar of type T (extracts T out of array<T>, if it's an array) */
	#define ISUNSPEC(T)	(ttrait<T>::isunspec)

	#if GPULOG_DEBUG && !__CUDACC__
	template<typename T>
	void dump_ttraits(const T& obj)
	{
		std::cerr << "align        = " << ALIGNOF(T) << "\n";
		std::cerr << "element size = " << ESIZE(T) << "\n";
		std::cerr << "dimen        = " << DIMEN(T) << "\n";
		std::cerr << "sizeof       = " << SIZEOF(T) << "\n";
		std::cerr << "isarray      = " << ISARRAY(T) << "\n";
		std::cerr << "isunspec     = " << ISUNSPEC(T) << "\n";
	}
	#endif

	} // namespace internal
} // namespace gpulog

#endif // bits_gpulog_ttraits_h__
