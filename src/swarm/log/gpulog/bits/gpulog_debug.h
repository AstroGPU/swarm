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

/*! \file gpulog_debug.h 
 *    \brief Defines debugging macros. 
 *
 *
 */

#ifndef bits_gpulog_debug_h__
#define bits_gpulog_debug_h__

//!
//! Debugging macros. Define GPULOG_DEBUG=1 to turn on
//!
#if GPULOG_DEBUG
	#define DBG(x) x
#else
	#define DBG(x)
#endif

#if !__CUDACC__
	#include <iostream>
	#include <cassert>
	#include <sstream>

	#define DHOST(x) DBG(x)
	#define DGPU(x)
	
	#define dev_assert(x) assert(x)
#else
	#define DHOST(x)	/* This hides the stuff that nvcc can't compile */

	#if __DEVICE_EMULATION__
		#define DGPU(x) DBG(x)
		// #define dev_assert(x) assert(x) /* has to be disabled on CUDA 2.3 or compilation fails */
		#define dev_assert(x) { if(!(x)) { printf("Assertion failed: " #x "\n"); exit(-1); } } /* This is a CUDA 2.3 compatible replacement for assert(); */
		// #define dev_assert(x)
	#else
		#define DGPU(x)
		#define dev_assert(x)
	#endif
#endif

#endif // bits_gpulog_debug_h
