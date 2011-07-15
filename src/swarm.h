/*************************************************************************
 * Copyright (C) 2010 by Mario Juric  and the Swarm-NG Development Team  *
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

/*! \file swarm.h
 *   \brief Main header file for swarm library
 *
 *   Declares class swarm_error.  Methods declared here are defined in swarmlib.cpp
 *
 *   Declares swarm_error, ensemble, cpu_ensemble, gpu_ensemble, writer (abstract), 
 *   Also contains some memory management funcctions
*/
#ifndef __swarm_h
#define __swarm_h

#include <stdexcept>
#include <string>
#include <cstring>
#include <map>
#include <cassert>
#include <cmath>
#include <vector>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cux/cux.h>

#include "ensemble.hpp"
#include "integrator.hpp"
#include "writer.h"
#include "stopwatch.h"
//#include "stopper.h"
//#include "propagator.h"

#ifdef THROW_IS_ABORT
	#include <cassert>
	#include <cstring>
        #include <cstdio>
#endif

#ifndef __CUDACC__ // CUDA 2.2 C++ bug workaround
	#include <sstream>
        #include <valarray>
#endif

#define CUDADEVICETOUSE 1

/// The main namespace for the Swarm-NG library
namespace swarm {

//class integrator;
//class writer;
//class gpu_ensemble;

typedef std::map<std::string, std::string> config;
typedef double real_time;
typedef float  real_mass;
typedef double real_pos;
typedef double real_vel;
typedef unsigned int uint;

/******** METHODS DEFINED IN swarmlib.cpp ***********/

/**
  \brief Configur grid for nthreads independent threads

   configure grid for nthreads independent threads, each requiring dynShmemPerThread of shared memory, and
   with each block needing staticShmemPerBlock of shared memory (usually to pass kernel invocation args.)
   The default for staticShmemPerBlock is reasonable (for small kernels), but not necessarily optimal.

   @param[out] gridDim
   @param[in] threadsPerBlock
   @param[in] nthreads
   @param[in] synShmemPerThread
   @param[in] staticShmemPerBlock
*/
bool configure_grid(dim3 &gridDim, int threadsPerBlock, int nthreads, int dynShmemPerThread = 0, int staticShmemPerBlock = 128);

/* Just prints message to cerr */
void debugger_stop();

/* Find the dimensions (bx,by) of a 2D grid of blocks that has as close
 * to nblocks blocks as possible*/
void find_best_factorization(unsigned int &bx, unsigned int &by, int nblocks);

/* Initialize the swarm library.  This function must be called before
 * any other */
void init(const config &cfg);

/// Load configuration from file fn
void load_config(config &cfg, const std::string &fn);

/// Load ensemble residing in files "name.XXX" where XXX \elem [0,nsys)
void load_ensemble(const std::string &name, cpu_ensemble &ens);


/******** METHODS DEFINED INLINE *************/

/*!
  \brief get a configuration value for 'key', throwing an error if it doesn't exist

  NOTE: heavy (unoptimized) function, use sparingly
*/
template<typename T>
void get_config(T &val, const config &cfg, const std::string &key)
{
	if(!cfg.count(key)) { ERROR("Configuration key '" + key + "' missing."); }
#ifndef __CUDACC__ // CUDA 2.2 C++ bug workaround
	std::istringstream ss(cfg.at(key));
	ss >> val;
#endif
}

/**
  \brief Typesafe re-allocator (convenience)
*/
template<typename T>
T* hostAlloc(T* var, int nelem, bool usePinned = false)
{
        if(!usePinned)
        {
                T* tmp = (T*)realloc(var, nelem*sizeof(T));
                if(tmp == NULL) ERROR("Out of host memory.");
                return tmp;
        }
        else
        {
                cudaThreadSynchronize();   // To prevent getting confused over other errors
                if(var!=NULL) hostFree(var);
                cudaError_t cudaMemStatus = cudaMallocHost((void**)&var,nelem*sizeof(T));
                if(cudaMemStatus!=cudaSuccess) ERROR(cudaGetErrorString(cudaMemStatus));
                return var;
        }
}

/**
  \brief Typesafe de-allocator (convenience)
*/
template<typename T>
void hostFree(T* var, bool usePinned = false)
{
        if(!usePinned)
        {
                ::free(var);
        }
        else
        {
                cudaThreadSynchronize();         // To prevent getting confused over other errors
                cudaError_t cudaMemStatus = cudaFreeHost(var);
                if(cudaMemStatus!=cudaSuccess) ERROR(cudaGetErrorString(cudaMemStatus));
        }
}

/*!
  \brief copy memory from host(CPU) to device(GPU)

  @param[out] dest destination to device momory 
  @param[in] src source from host memory
*/
template<typename T>
inline void memcpyToGPU(T *dest, const T *src, int nelem = 1)
{
	cuxErrCheck( cudaMemcpy(dest, src, nelem*sizeof(T), cudaMemcpyHostToDevice) );
}

/*!
  \brief copy memory from device(GPU) to host(CPU)

  @param[out] dest destination to host momory 
  @param[in] src source from device memory
*/
template<typename T>
inline void memcpyToHost(T *dest, const T *src, int nelem = 1)
{
	cuxErrCheck( cudaMemcpy(dest, src, nelem*sizeof(T), cudaMemcpyDeviceToHost) );
}

/********** UTILITIES *************/

// NOTE: The ifdef here is a workaround for CUDA 2.2 device emulation mode bug, where C++
// is disabled in emulation mode. If you want to use the functions below, use them only
// from units compiled with g++
#ifndef __CUDACC__
	void calc_total_energy(const cpu_ensemble &ens, std::valarray<double> &E);
#endif

} // end namespace swarm

#include <iostream>
#define $$(x) (std::cerr << __FILE__ << "(" << __FUNCTION__ << "):" << __LINE__ << " |> " << (x) << std::endl)
#define $_(x) (std::cerr << __FILE__ << "(" << __FUNCTION__ << "):" << __LINE__ <<  " " << (#x) << " = " << (x) << std::endl)
#define $$$ (std::cerr << __FILE__ << "(" << __FUNCTION__ << "):" << __LINE__ << " @@ " << std::endl)
#endif
