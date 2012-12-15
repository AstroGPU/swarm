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

/*! \file gpulog_ilogstream.h 
 *    \brief Handles log stream records. 
 *
 *
 */

#ifndef bits_gpulog_ilogstream_h__
#define bits_gpulog_ilogstream_h__

namespace gpulog
{
namespace internal
{

	struct unexpected_end_of_buffer_exception : std::exception {
		const char * what() {
			return "Unexpected end of buffer, the record retrieved is only half correct";
		}
	};

        //! a stream of logrecords
	struct ilogstream
	{
	protected:
		header hend;

	protected:
		const char *const ptr;
		size_t at, len;

	public:
	        //! Class constructors for different stream input
		__device__ __host__ ilogstream(const char *ptr_, size_t len_) : ptr(ptr_), at(0), len(len_), hend(-1, 0) {}
	        //! Class constructors for different stream input
		__device__ __host__ ilogstream(const host_log &l) : ptr(l.internal_buffer()), at(0), len(l.size()), hend(-1, 0) {}

	        //! Class constructors for different stream input
		__device__ __host__ logrecord next()
		{
			// if we're at the end
			if(at == len) { return logrecord((char *)&hend); }

			// return next packet
			logrecord rec(ptr + at);
			at += rec.len();

			if(at > len)
				throw unexpected_end_of_buffer_exception();

			return rec;
		}

		const size_t& cursor(){
			return at;
		}

	};

}
}

#endif
