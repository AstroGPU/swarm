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

//#include "config.h"

#include <iostream>
#include <fstream>

//#include <astro/exceptions.h>
#include <astro/util.h>
#include <astro/macros.h>
//#include <astro/math.h>
//#include <astro/system/log.h>
//#include <astro/io/format.h>
//#include <astro/useall.h>

#include <vector>
//#include "gpu.h"
//#include "io.h"
//#include "analysis.h"
//#include "spline.h"

#include <cuda.h>
#include <cuda_runtime.h>

#include "cux.h"
#include "cuda_rng.h"

#define ASSERT assert
#define MLOG(x) std::cerr
#define DLOG(x) std::cerr

size_t arrayMemSize_impl(size_t nx, size_t ny, size_t nz, size_t align, size_t elementSize)
{
	size_t size;	// size of the area to allocate, in bytes

	if(ny == 1 && nz == 1)
	{
		// no extra alignment padding for 1D array
		size = elementSize*nx;
	}
	else
	{
		size = roundUpModulo(nx*elementSize, align) * ny * nz;
	}

	return size;
}

cuxSmartPtr_impl_t::cuxSmartPtr_impl_t(size_t es, size_t pitch, size_t width, size_t height, size_t depth)
{
	ASSERT(pitch >= width*es);
	ASSERT(depth >= 1);
	ASSERT(height >= 1);

	// array layout
	m_data.extent[0] = pitch;
	m_data.extent[1] = height;
	m_data.extent[2] = depth;
	m_width = width;
	m_elementSize = es;

	// reference counting and garbage collection
	refcnt = 1;
	all_cuxSmartPtrs.insert(this);

	// NOTE: storage is lazily allocated the first time this pointer
	// is accessed through syncTo* methods
	m_data.ptr = NULL;
	slave = NULL;
	cuArray = NULL;
	onDevice = false;
	cleanCudaArray = false;
}

cuxSmartPtr_impl_t::~cuxSmartPtr_impl_t()
{
	ASSERT(boundTextures.empty()); // Assert the textures were unbound, before gc() runs

	// delete all slave copies
	gc();

	// delete the master
	if(!onDevice)
	{
		delete [] m_data.ptr;
	}
	else
	{
		cuxErrCheck( cudaFree(m_data.ptr) );
	}

	all_cuxSmartPtrs.erase(this);
}

cuxSmartPtr_impl_t::allocated_pointers::~allocated_pointers()
{
	if(!empty())
	{
		MLOG(verb1) << "ERROR: Memory leak -- " << size() << " cuxSmartPtr<> pointers were not deallocated\n";
		int k = 0;
		FOREACH(*this)
		{
			cuxSmartPtr_impl_t *p = *i;
			MLOG(verb1) << "  " << k << ": " << p->m_width << " x " << p->m_data.extent[1] << " x " << p->m_data.extent[2]
				<< " (elementSize=" << p->m_elementSize << ")";
			k++;
		}
	}
}

cuxSmartPtr_impl_t::allocated_pointers cuxSmartPtr_impl_t::all_cuxSmartPtrs;
void cuxSmartPtr_impl_t::global_gc()
{
	FOREACH(all_cuxSmartPtrs)
	{
		(*i)->gc();
	}
}

void cuxSmartPtr_impl_t::gc()
{
	// delete slave (global memory) copy
	if(onDevice)
	{
		delete [] slave;
	}
	else
	{
		cuxErrCheck( cudaFree(slave) );
	}
	slave = NULL;

	// if the cudaArray is dirty, or there are no textures bound to it
	// assume it's available for deletion
	if(!cleanCudaArray || boundTextures.empty())
	{
		if(cuArray)
		{
			cuxErrCheck( cudaFreeArray(cuArray) );
		}
		cuArray = NULL;
		cleanCudaArray = false;
	}
}

void *cuxSmartPtr_impl_t::syncTo(bool device)
{
	if(onDevice != device)
	{
		std::swap(slave, m_data.ptr);
	}

	// Allocate m_data.ptr if needed.
	if(!m_data.ptr)
	{
		if(device)
		{
			//cuxErrCheck( cudaMalloc((void**)&m_data.ptr, memsize()) );
			cudaError err = cudaMalloc((void**)&m_data.ptr, memsize());
			if(err == cudaErrorMemoryAllocation)
			{
				global_gc();
				err = cudaMalloc((void**)&m_data.ptr, memsize());
			}
			cuxErrCheck(err);
		}
		else
		{
			m_data.ptr = new char[memsize()];
			//memset(m_data.ptr, 0xff, memsize());	// debugging
		}
	}

	// copy slave -> m_data.ptr (if there's something to copy)
	if(onDevice != device && slave)
	{
		cudaMemcpyKind dir = device ? cudaMemcpyHostToDevice : cudaMemcpyDeviceToHost;
		cuxErrCheck( cudaMemcpy(m_data.ptr, slave, memsize(), dir) );
	}

	onDevice = device;

	// assume the sync dirtied up the textures
	cleanCudaArray = false;

//	gc(); // agressive garbage collection while debugging

	return m_data.ptr;
}

#define GC_AND_RETRY_IF_FAIL(x) \
	{ \
		cudaError err = (x); \
		if(err == cudaErrorMemoryAllocation) \
		{ \
			global_gc(); \
			err = (x); \
		} \
		cuxErrCheck(err); \
	}

cudaArray *cuxSmartPtr_impl_t::getCUDAArray(cudaChannelFormatDesc &channelDesc)
{
	ASSERT(channelDesc.x + channelDesc.y + channelDesc.z + channelDesc.w == m_elementSize*8);

	// FIXME: This all seems to be majorly fu*ked up, as CUDA devemu
	// has bugs with cudaMalloc3DArray that has any of the extent dimensions
	// set to zero. Will have to be fixed by trial-and-error on the real GPU.
	if(!cleanCudaArray)
	{
		syncToHost();	// ensure the data is on the host

		// array size in elements (we're going to need this later)
		cudaExtent ex = make_cudaExtent(m_width, m_data.extent[1], m_data.extent[2]);
		ASSERT(ex.width > 0);

		if(!cuArray)	// allocate if needed
		{
			if(ex.depth > 1)
			{
				// 3D arrays
				GC_AND_RETRY_IF_FAIL( cudaMalloc3DArray(&cuArray, &channelDesc, ex) );
			}
			else
			{
				// 2D and 1D arrays
				GC_AND_RETRY_IF_FAIL( cudaMallocArray(&cuArray, &channelDesc, ex.width, ex.height) );
			}
		}

		// copy
		if(ex.depth > 1)
		{
			// 3D arrays
			cudaMemcpy3DParms par = { 0 };
			par.srcPtr = make_cudaPitchedPtr(m_data.ptr, m_data.extent[0], ex.width, ex.height);
			par.dstArray = cuArray;
			par.extent = ex;
			par.kind = cudaMemcpyHostToDevice;
			cuxErrCheck( cudaMemcpy3D(&par) );
		}
		else
		{
			// 2D and 1D arrays
			cuxErrCheck( cudaMemcpy2DToArray(cuArray, 0, 0, m_data.ptr, m_data.extent[0], ex.width*m_elementSize, ex.height, cudaMemcpyHostToDevice) );
		}

		cleanCudaArray = true;
	}

	ASSERT(cuArray);
	return cuArray;
}

// texture access
void cuxSmartPtr_impl_t::bind_texture(textureReference &texref)
{
	cudaArray *cuArray = getCUDAArray(texref.channelDesc);
	cuxErrCheck( cudaBindTextureToArray(&texref, cuArray, &texref.channelDesc) );
	
	boundTextures.insert(&texref);
}

void cuxSmartPtr_impl_t::unbind_texture(textureReference &texref)
{
	cuxErrCheck( cudaUnbindTexture(&texref) );

	ASSERT(boundTextures.count(&texref));
	boundTextures.erase(&texref);
}

//
// texture load and creation utilities
//

cuxTexture<float, 1> load_constant_texture_1D(float val, float X0, float X1)
{
	cuxTexture<float, 1> tex(2);
	tex(0U) = val;
	tex(1U) = val;
	tex.coords[0].x = X0;
	tex.coords[0].y = 1./(X1 - X0);
	return tex;
}

cuxTexture<float, 3> load_constant_texture_3D(
	float val,
	float x0, float x1,
	float y0, float y1,
	float z0, float z1
)
{
	float2 tcx = make_float2(x0, 1./(x1-x0)), tcy = make_float2(y0, 1./(y1-y0)), tcz = make_float2(z0, 1./(z1-z0));
	cuxTexture<float, 3> tex(2, 2, 2, tcx, tcy, tcz);

	FORj(i, 0, 2) FORj(j, 0, 2) FORj(k, 0, 2)
		tex(i, j, k) = val;

	return tex;
}

//
// Return texcoord that will map x to imgx and y to imgy
//
float2 texcoord_from_range(float imgx, float imgy, float x, float y)
{
	float2 tc;

	tc.x = (-imgy * x + imgx * y)/(imgx - imgy);
	tc.y = (imgx - imgy) / (x - y);

	return tc;
}

/***********************************************************************/

// CUDA emulation for the CPU
// Used by CPU versions of CUDA kernels
__TLS char impl_shmem[16384];
namespace gpuemu	// prevent collision with nvcc's symbols
{
	__TLS uint3 blockIdx;
	__TLS uint3 threadIdx;
	__TLS uint3 blockDim;	// Note: uint3 instead of dim3, because __TLS variables have to be PODs
	__TLS uint3 gridDim;		// Note: uint3 instead of dim3, because __TLS variables have to be PODs
}

__TLS int  active_compute_device;


///////////////////////////////////////////////////////////
// CUDA helpers

const char *cpuinfo()
{
	static char buf[1000];
	FILE *f = popen("cat /proc/cpuinfo | grep 'model name' | head -n 1 | awk -F': ' '{ print $2}'", "r");
	fgets(buf, 1000, f);
	pclose(f);

	int len = strlen(buf);
	if(len && buf[len-1] == '\n') buf[len-1] = 0;
	return buf;
}


void abort_on_cuda_error(cudaError err)
{
	if(err == cudaSuccess) { return; }

	MLOG(verb1) << "CUDA ERROR: " << cudaGetErrorString(err);
	//abort();
	exit(-100);
}

void cuxErrCheck_impl(cudaError err, const char *fun, const char *file, const int line)
{
	if(err != cudaSuccess)
	{
		MLOG(verb1) << "CUDA ERROR: In " << fun << " (" << file << ":" << line << ")\n";
		abort_on_cuda_error(err);
//		throw cuxException(err);
	}
}

static int cuda_initialized = 0;	// CUDA initialized
static int cuda_enabled = 0;		// Should GPU acceleration be used?

/**
	Initialize cux library and CUDA device.
*/
bool cux_init()
{
	if(cuda_initialized) { return true; }

	// get requested device from environment
	int dev;
	const char *devStr = getenv("CUDA_DEVICE");
	bool autoselect = devStr == NULL;

	if(!autoselect)
	{
		dev = atoi(devStr);

		// disable GPU acceleration
		if(dev == -1)
		{
			cuda_initialized = 1;
			cuda_enabled = 0;

			MLOG(verb1) << "GPU accelerator: Using CPU: \"" << cpuinfo() << "\"";
			return true;
		}
		else
		{
			cuxErrCheck( cudaSetDevice(dev) );
		}
	}

#if !CUDA_DEVEMU
	// ensure a CUDA context is created and fetch the active
	// device id
	void *tmp;
	cuxErrCheck( cudaMalloc(&tmp, 1024) );
	cuxErrCheck( cudaFree(tmp) );
	cuxErrCheck( cudaGetDevice(&dev) );
#endif

#if !CUDA_DEVEMU
	// get device properties
	cudaDeviceProp deviceProp;
	cuxErrCheck( cudaGetDeviceProperties(&deviceProp, dev) );

	char buf[1000];
	snprintf(buf, sizeof(buf), "GPU accelerator: Using Device %d: \"%s\"%s", dev, deviceProp.name, autoselect ? " (autoselected)" : "");
	MLOG(verb1) << buf;
#else
	MLOG(verb1) << "GPU accelerator: Using Device Emulation";
#endif

#if !CUDA_DEVEMU
	// Memory info
	unsigned free = 0, total = 0;
	cuxErrCheck( (cudaError)cuMemGetInfo(&free, &total) );
	MLOG(verb2) << "Device memory (free, total): " << free / (1<<20) << "M, " << total / (1<<20) << "M" << "\n";
#endif

	cuda_initialized = 1;
	cuda_enabled = 1;

	return true;
}
