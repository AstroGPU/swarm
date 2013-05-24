/*************************************************************************
 * Copyright (C) 2011 by Saleh Dindar and the Swarm-NG Development Team  *
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

/*! \file common.hpp
 *   \brief Common external (and some internal) library headers used by all components of swarm.
 *
 *   This file can be compiled as a precompiled header.
 *
 *   External libraries referenced:
 *    - Standard C++ Library
 *    - Boost libraries (bind, shared_ptr)
 *    - CUDA runtime library
 *    - Linux libc headers
 *
 *   Internal headers referenced:
 *    - Common error handling: base class for all other exceptions
 *
 *   Some useful debugging and CUDA macros
 *
 *
*/
#pragma once

#ifdef __CUDACC__
	// GCC 4.7 bug workaround
	#undef _GLIBCXX_ATOMIC_BUILTINS
	#undef _GLIBCXX_USE_INT128
#endif


// Standard C++ Library
#include <cassert>
#include <cmath>
#include <stdint.h>
#include <cstring>
#include <cstdio>

// STL
#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <map>
#include <vector>
#include <set>
#include <algorithm>
#include <limits>
#ifndef __CUDACC__ 

	// CUDA 2.2 C++ bug workaround
	#include <sstream>
	#include <valarray>
#endif

#ifndef __CUDACC__
// Boost Libraries
#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>

using boost::shared_ptr;

// CUDA libraries
#include <cuda.h>
#include <cuda_runtime.h>
#endif


// POSIX headers
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>

#include "runtime_error.hpp"

//
// Struct alignment is handled differently between the CUDA compiler and other
// compilers (e.g. GCC, MS Visual C++ .NET)
//
#ifdef __CUDACC__
	#define ALIGN(x)  __align__(x)
#else
	#if defined(_MSC_VER) && (_MSC_VER >= 1300)
		// Visual C++ .NET and later
		#define ALIGN(x) __declspec(align(x))
	#else
		#if defined(__GNUC__)
			// GCC
			#define ALIGN(x)  __attribute__ ((aligned (x)))
		#else
		// all other compilers
			#define ALIGN(x)
		#endif
	#endif
#endif

//! To debug output a value
#define $$(x) (std::cerr << __FILE__ << "(" << __FUNCTION__ << "):" << __LINE__ << " |> " << (x) << std::endl)

#define $PRINT(x) (std::cerr << __FILE__ <<  __LINE__ << ": " << x << std::endl)

//! Show the an expression and its value for debugging
#define $_(x) (std::cerr << __FILE__ << "(" << __FUNCTION__ << "):" << __LINE__ <<  " " << (#x) << " = " << (x) << std::endl)

//! To indicate when the execution reaches a specific line in code
#define $$$ (std::cerr << __FILE__ << "(" << __FUNCTION__ << "):" << __LINE__ << " @@ " << std::endl)
#define SYNC cudaThreadSynchronize()

