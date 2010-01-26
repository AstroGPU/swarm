/***************************************************************************
 *   Copyright (C) 2004 by Mario Juric                                     *
 *   mjuric@astro.Princeton.EDU                                            *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
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

#ifndef cux_lowlevel__
#define cux_lowlevel__

// mjuric: short-circuit CUDA auto-detection
#define HAVE_CUDA 1

//////////////////////////////////////////////////////////////////////////
// CUDA (or emulation) APIs
//////////////////////////////////////////////////////////////////////////
#if HAVE_CUDA
	#include <cuda_runtime.h>
#else
	#include "cuda_emulation.h"
#endif

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

// Thread local storage -- use for shared memory emulation in CPU mode
//#define __TLS __thread
#define __TLS

//////////////////////////////////////////////////////////////////////////
// Shared memory access and CPU emulation
//////////////////////////////////////////////////////////////////////////
#if __CUDACC__
	extern __shared__ char impl_shmem[];
#else
	// For CPU versions of GPU algorithms
	extern __TLS char impl_shmem[16384];
	namespace gpuemu // prevent collision with nvcc's symbols
	{
		extern __TLS uint3 blockIdx;
		extern __TLS uint3 threadIdx;
		extern __TLS uint3 blockDim;		// Note: uint3 instead of dim3, because __TLS variables have to be PODs
		extern __TLS uint3 gridDim;		// Note: uint3 instead of dim3, because __TLS variables have to be PODs
	}
	using namespace gpuemu;
#endif

#define shmem(type) ((type*)impl_shmem)

bool calculate_grid_parameters(dim3 &gridDim, int threadsPerBlock, int neededthreads, int dynShmemPerThread, int staticShmemPerBlock);

#endif
