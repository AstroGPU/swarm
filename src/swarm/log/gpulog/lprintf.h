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

#ifndef lprintf_h__
#define lprintf_h__

#include <cstring>
#include "gpulog.h"

	/*
		printf API
	*/
#if 0
	/* CUDA 2.3 compatible float->double converter */
	template<typename T> __device__ __host__ inline  const T& f2d_fun(const T &v)      { return v; }
	__device__ inline double f2d_fun(const float &v)		// specialization for cast of float to double (used by printf)
	{
		return v;
	}
	
	#define f2d(T, x) f2d_fun(x)
#else
	/* CUDA 2.2 compatible float->double converter hack */
	template<typename T> struct f2dstruct
	{
		typedef const T& cTref;
		__host__ __device__ static cTref cvt(cTref v) { return v; }
	};
	template<> struct f2dstruct<float>
	{
		__host__ __device__ static double cvt(const float &v) { return v; }
	};

	#define f2d(T, x) f2dstruct<T>::cvt(x)
#endif

	#if 0
	template<typename L, int N, typename T1, typename T2, typename T3>
		__host__ __device__ inline void lprintf(L &log, const char (&fmt)[N], const T1 &v1, const T2 &v2, const T3 &v3)
		{
			log.write(MSG_PRINTF, fmt, f2d(v1), f2d(v2), f2d(v3));
		}
	#else
		#include "bits/gpulog_printf.h"
	#endif


	#if !__CUDACC__
	inline std::string run_printf(gpulog::logrecord &lr)
	{
        	using namespace gpulog;

		std::string res;
		if(lr.msgid() != MSG_PRINTF)
		{
			return "";
		}

		// slurp up the format string
		char fmtbuf[1024], *fmt = fmtbuf;
		lr >> fmt;

		double data_double;
		int data_int;
		char data_str[1024];

		// Now run through it, printing everything we can. We must
		// run to every % character, extract only that, and use printf
		// to format it.
		char buf[1024];
		std::ostringstream out;
		char *p = strchr ( fmt, '%' );
		while ( p != NULL )
		{
			// Print up to the % character
			*p = '\0';
			out << fmt;
			*p = '%';           // Put back the %

			// Now handle the format specifier
			char *format = p++;         // Points to the '%'
			p += strcspn ( p, "%cdiouxXeEfgGaAnps" );
			if ( *p == '\0' )           // If no format specifier, print the whole thing
			{
				fmt = format;
				break;
			}

			char specifier = *p++;
			char c = *p;        // Store for later
			*p = '\0';
			switch ( specifier )
			{
					// These all take integer arguments
				case 'c':
				case 'd':
				case 'i':
				case 'o':
				case 'u':
				case 'x':
				case 'X':
				case 'p':
					lr >> data_int;
					sprintf(buf, format, data_int);
					out << buf;
					break;

					// These all take double arguments
				case 'e':
				case 'E':
				case 'f':
				case 'g':
				case 'G':
				case 'a':
				case 'A':
					lr >> data_double;
					sprintf(buf, format, data_double);
					out << buf;
					break;

					// Strings are handled in a special way
				case 's':
					lr >> data_str;
					sprintf(buf, format, data_str);
					out << buf;
					break;

					// % is special
				case '%':
					out << "%%";
					break;

					// Everything else is just printed out as-is
				default:
					out << format;
					break;
			}
			*p = c;                     // Restore what we removed
			fmt = p;                    // Adjust fmt string to be past the specifier
			p = strchr(fmt,'%');    // and get the next specifier
		}

		// Print out the last of the string
		out << fmt;
		return out.str();
	}

	inline int replay_printf(std::ostream &out, gpulog::ilogstream &ls)
	{
		using namespace gpulog;

		logrecord lr;
		int count = 0;
		while(lr = ls.next())
		{
			if(lr.msgid() != MSG_PRINTF) { continue; }

			// process the printf...
			out << run_printf(lr);
			count++;
		}
		return count;
	}

	inline int replay_printf(std::ostream &out, const gpulog::host_log &log)
	{
		gpulog::ilogstream ls(log);
		return replay_printf(out, ls);
	}
	#endif



#endif // lprintf_h__
