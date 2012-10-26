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

/*! \file gpulog_types.h 
 *    \brief Defines the interface for logging data types. 
 *
 *
 */

#ifndef bits_gpulog_types_h__
#define bits_gpulog_types_h__

#include <cuda_runtime.h>

namespace gpulog
{
	namespace internal
	{
		// unbound array alloc request class
		template<typename T>
		struct array
		{
			int nelem;	// requested size (in elements)

			__device__ inline array(int n_) : nelem(n_) {}
		};

		// arginfo (mainly used for debugging)
		struct arginfo
		{
			int arg;			// which arg is it
			int align, size, dim;		// scalar type alignment, scalar type size, dimension (if presized array)
			bool isarray; int nelem;	// is array<>, number of allocated elements
			int begin, len;			// beginning offset and total length (bytes, w/o padding)
		};

		// log record header
		struct header
		{
			int msgid;	// record ID of the log record
			int len;	// byte length of the data (including all padding and the size of the header)

			#if ARGINFO
			int nargs;	// number of arguments (including the header)
			int infos;	// offset to array of arginfo structures
			#endif

			__host__ __device__ int dataoffset() const
			{
				int offs = sizeof(*this);
				#if ARGINFO
				offs += nargs * sizeof(arginfo);
				#endif
				return offs;
			}

			__host__ __device__ static inline void new_header(header &h, int msgid_, int len_)
			{
				h.msgid = msgid_; h.len = len_;
				#if ARGINFO
				h.nargs = 0;
				h.infos = 0;
				#endif
			}

			__host__ __device__ inline header(int msgid_, int len_) { new_header(*this, msgid_, len_); }
		};

		// "unspecified datatype" marker for pktsize structure
		struct Tunspec { };

		// datatype with maximum allowed alignment
		struct ALIGN(16) Tmaxalign {};
	}

}

#endif // bits_gpulog_types_h__
