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

#ifndef cuda_cux__
#define cuda_cux__

#include <assert.h>
#include <stdint.h>
#include <map>
#include <set>
#include <algorithm>

#include "cux_lowlevel.h"

/*
	Initialize compute device. Must be called before any other
	calls involving the device.
*/
bool cux_init();

//////////////////////////////////////////////////////////////////////////
// Useful host and device functions
//////////////////////////////////////////////////////////////////////////

/*
	Rounds up v to the nearest integer divisible by mod. Usually used to
	compute size of padded and aligned arrays.
*/
inline int roundUpModulo(int v, int mod)
{
	int r = v % mod;
	int pitch = r ? v + (mod-r) : v;
	return pitch;
}

/**
	Computes the global linear ID of the thread. Used from kernels.
*/
inline __device__ uint32_t threadID()
{
	// NOTE: Supports 3D grids with 1D blocks of threads
	const uint32_t id = ((blockIdx.z * gridDim.y + blockIdx.y) * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;
	return id;
}

/** ********************************************************************
    cux API classes
************************************************************************/

/**
 * cuxException
 * \brief class to be thrown if a CUDA-related error is detected
*/
struct cuxException
{
	cudaError err;
	cuxException(cudaError err_) : err(err_) {}
	const char *msg() const { return cudaGetErrorString(err); }
};

/*
	cuxErrCheck macro -- aborts with message if the enclosed call returns != cudaSuccess
*/
void cuxErrCheck_impl(cudaError err, const char *fun, const char *file, const int line);
#define cuxErrCheck(expr) \
	cuxErrCheck_impl(expr, __PRETTY_FUNCTION__, __FILE__, __LINE__)


/*
	Compute the size, in bytes, of an aligned (strided) 1/2/3D array.
*/
size_t arrayMemSize_impl(size_t nx, size_t ny, size_t nz, size_t align, size_t elementSize);
template<typename T>
	size_t arrayMemSize(size_t nx, size_t ny = 1, size_t nz = 1, size_t align = 128)
	{
		return arrayMemSize_impl(nx, ny, nz, align, sizeof(T));
	}

/*
	cuxNew<T> -- Allocate an array on the current device

	Low level API -- cuxDevicePtr<> should be used instead.

	If 2D or 3D array is requested, the first dimension will be byte aligned
	according to the align argument (default: 128 bytes).

	Hides the raw CUDA API.
*/
template<typename T>
	T* cuxNew(size_t nx, size_t ny = 1, size_t nz = 1, size_t align = 128)
	{
		size_t size = arrayMemSize<T>(nx, ny, nz, align);

		T *devptr;
		cuxErrCheck( cudaMalloc((void**)&devptr, size) );
		return devptr;
	}

/*
	cuxFree() -- Free memory on the current device

	Low level API -- cuxDevicePtr<> should be used instead.

	Hides the raw CUDA API.
*/
template<typename T>
	void cuxFree(T *v)
	{
		cuxErrCheck( cudaFree(v) );
	}


/*
	cuxUploadConst -- Copy a variable to __constant__ memory variable

	The __constant__ variable can be given via address (if compiling with nvcc)
	or a string name (if called from gcc-compiled code)

	Hides the raw CUDA API.
*/

template<typename T>
	void cuxUploadConst(T &dest, const T &source)
	{
		unsigned size = sizeof(source);
#ifdef __CUDACC__
		cuxErrCheck( cudaMemcpyToSymbol(dest, &source, size) );
#elif BUILD_FOR_CPU
		memcpy(&dest, &source, size);
#else
		assert(0);
		//#error cuxUploadConst can be used only in .cu files, or when BUILD_FOR_CPU is defined
#endif
	}

template<typename T>
	void cuxUploadConst(const char *symbol, const T &source)
	{
		unsigned size = sizeof(source);
		//fprintf(stderr, "cuxUploadConst: Uploading %u bytes to symbol %s.\n", size, symbol);
		cuxErrCheck( cudaMemcpyToSymbol(symbol, &source, size) );
	}

/*
	arrayPtr -- Multidimensional strided array interface.

	Size: sizeof(T*) + sizeof(size_t)*(dim-1)

	Stores the pointer to the array, and the size (in bytes) of first dim-1
	dimensions. If dim > 1, ensures the first dimension occupies is a
	multiple of aling bytes. If dim == 1, it becomes equivalent (in size
	and usage) to a simple C pointer.

	The user is responsible for memory allocation/deallocation, using the
	make_arrayPtr() functions (see below), or using cuxDevicePtr<> class (below).

	The array may be in either host or device memory and must be accessed
	accordingly; arrays in device memory can only be accessed in kernel
	code, arrays in host memory can only be accessed from host code.
*/
template<typename T, int dim = 1>
	struct arrayPtr
	{
		T *ptr;
		uint32_t extent[dim-1];	// extent[0] == width of a row of data, in bytes
					// extent[1] == number of rows of data in a slice of a 3D data cube
					// extent[2] == number of slices in a 3D data sub-cube of a 4D data cube (etc...)

		// Access
#if 0 // Oh, there are bigger problems about allocating memory, so this won't be such a trivial upgrade 
		__host__ __device__ T &operator()(const uint32_t w,  const uint32_t x, const uint32_t y, const uint32_t z) const	// 4D accessor 
		{
			return *((T*)((char*)ptr + x * extent[0] + y * extent[0] * extent[1] + z * extent[0] * extent[1] * extent[2]) + w);
		}
#endif
		__host__ __device__ T &operator()(const uint32_t x, const uint32_t y, const uint32_t z) const	// 3D accessor (valid only if dim >= 3)
		{
			return *((T*)((char*)ptr + y * extent[0] + z * extent[0] * extent[1]) + x);
		}
		__host__ __device__ T &operator()(const uint32_t x, const uint32_t y) const	// 2D accessor (valid only if dim >= 2)
		{
			return *((T*)((char*)ptr + y * extent[0]) + x);
		}
		__host__ __device__ T &operator()(const uint32_t i) const			// 1D accessor (valid only if dim >= 1)
		{
			return ptr[i];
		}
		__host__ __device__ T &operator[](const uint32_t i) const			// 1D accessor (valid only if dim >= 1)
		{
			return ptr[i];
		}
		__host__ __device__ operator bool() const					// comparison with NULL
		{
			return ptr != NULL;
		}
		__host__ __device__ void reset()						// allow the setting of pointer to NULL
		{
			ptr = NULL;
		}
	};

/*
	cuxDevicePtr<T, dim, align> -- n-dimensional array on the device

	Size: sizeof(T*) + sizeof(size_t)*(dim-1)

	Derived from arrayPtr. Encapsulates a pointer to array on the device,
	simplifying upload, download, allocation and deallocation and abstracting
	it from the CUDA APIs.

	Provides NO automatic memory management (e.g., RAII); must be allocated
	and deallocated manually via member functions. For a smart pointer, see
	cuxSmartPtr (below).
*/

template<typename T, int dim = 1, int align = 128>
	struct cuxDevicePtr : public arrayPtr<T, dim>
	{
	protected:
		size_t memsize(size_t lastdim)
		{
			size_t size = lastdim;
			if(dim == 1) return size*sizeof(T);

			for(int i = 0; i != dim-1; i++)
			{
				size *= this->extent[i];
			}
			return size;
		}
	public:
		void upload(const T* src, size_t lastdim)
		{
			if(!this->ptr)
			{
				if(dim == 1) { alloc(lastdim); }	// auto-allocation allowed only for 1D arrays (convenience)
				else { assert(this->ptr); }
			}

			size_t size = memsize(lastdim);
			cuxErrCheck( cudaMemcpy(this->ptr, src, size, cudaMemcpyHostToDevice) );
		}

		void download(T* dest, size_t lastdim)
		{
			if(!this->ptr)
			{
				if(dim == 1) { alloc(lastdim); }	// auto-allocation allowed only for 1D arrays (convenience)
				else { assert(this->ptr); }
			}

			size_t size = memsize(lastdim);
			cuxErrCheck( cudaMemcpy(dest, this->ptr, size, cudaMemcpyDeviceToHost) );
		}

		void realloc(size_t nx, size_t ny = 1, size_t nz = 1)
		{
			free();
			alloc(nx, ny, nz);
		}
		void alloc(size_t nx, size_t ny = 1, size_t nz = 1)
		{
			assert(dim <= 4); // not implemented for dim > 4

			if(dim >  1) { this->extent[0] = roundUpModulo(nx*sizeof(T), align); }
			if(dim >  2) { this->extent[1] = ny; }
			if(dim >  3) { this->extent[2] = nz; }

			this->ptr = cuxNew<T>(nx, ny, nz, align);
		}
		void free()
		{
			cuxFree(this->ptr);
			this->ptr = NULL;
		}

		cuxDevicePtr<T, dim, align> operator =(T *p) { this->ptr = p; return *this; }
		cuxDevicePtr<T, dim, align> operator =(const cuxDevicePtr<T> &p) { this->ptr = p.ptr; return *this; }
	};

/*
	cuxDeviceAutoPtr<T, dim, align> -- n-dimensional array on the device
		with automatic deallocation in destructor
**/
template<typename T, int dim = 1, int align = 128>
	struct cuxDeviceAutoPtr : public cuxDevicePtr<T, dim, align>
	{
		cuxDeviceAutoPtr()
		{
			this->ptr = NULL;
		}
		cuxDeviceAutoPtr(int nx, int ny = 1, int nz = 1)
		{
			this->alloc(nx, ny, nz);
		}
		~cuxDeviceAutoPtr()
		{
			this->free();
		}
		operator T*() const { return this->ptr; }
		void get(T* val, int n = 1)
		{
			assert(dim == 1); // TODO: implement for higher-D arrays
			memcpyToHost(val, this->ptr, n);
		}
		void memset(int val, int n = 1)
		{
			assert(dim == 1); // TODO: implement for higher-D arrays
			cudaMemset(this->ptr, val, sizeof(T)*n);
		}
	};


////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

// Pointer to an n-D array in device memory, obtained from cuxSmartPtr
//	- Used to obtain a pointer to and access device memory (in kernels)
//	- Should be passed to device kernels
//	- Thin veneer over arrayPtr, to facilitate obtaining by cast from cuxSmartPtr<>
//	- Must be obtained via cast from cuxSmartPtr<>
//	- Remains valid until a hptr<> or a texture_bind is called for the parent cuxSmartPtr<>
//	- Must not be manually deallocated
template<typename T, int dim = 2>
struct gptr : public arrayPtr<T, dim>
{
};

// Pointer to n-D array in host memory, obtained from cuxSmartPtr
//	- Used to obtain a pointer to and access host memory
//	- Thin veneer over arrayPtr, to facilitate obtaining by cast from cuxSmartPtr<>
//	- Must be obtained via cast from cuxSmartPtr<>
//	- Remains valid until a gptr<> or a texture_bind is called for the parent cuxSmartPtr<>
//	- Must not be manually deallocated
template<typename T, int dim = 2>
struct hptr : public arrayPtr<T, dim>
{
};

// Inner part, low level implementation of cuxSmartPtr<T>. See the documentation of cuxSmartPtr<T>
// for more details.
//     - only meant to be used from cuxSmartPtr
//     - points to a well-formed block 3D block of elements of size m_elementSize,
//	 with byte dimensions in m_data.extent[0-2], and logical width given in
//	 m_width.
//     - is reference counted, auto de-allocates on all devices upon final release
struct cuxSmartPtr_impl_t
{
public:
	// data members
	arrayPtr<char, 4> m_data;	// master copy of the data (can be on host or device, depending on onDevice)
	char *slave;			// slave copy of the data  (can be on host or device, depending on onDevice)
	cudaArray* cuArray;		// CUDA array copy of the data

	bool onDevice;			// true if the "master copy" of the data is on the device
	bool cleanCudaArray;		// true if the last access operation was obtaining a reference to cudaArray

	uint32_t m_elementSize;		// size of the array element (in bytes)
	uint32_t m_width;		// logical width of the array (in elements)

	int refcnt;			// number of cuxSmartPtrs pointing to this m_impl

	std::set<textureReference *> boundTextures;	// list of textures bound to the cuArray

public:
	// constructors/destructors
	cuxSmartPtr_impl_t(size_t elementSize, size_t pitch, size_t width, size_t height = 1, size_t depth = 1);
	~cuxSmartPtr_impl_t();

	// aux methods
	operator bool() const					// true if the pointer is considered non-null (effectively defines what non-null means)
	{
		return m_data.extent[0] != 0;
	}
	uint32_t memsize() const				// number of bytes allocated
	{
		uint32_t size = m_data.extent[0]*m_data.extent[1]*m_data.extent[2];
		return size;
	}

	// reference management
	cuxSmartPtr_impl_t *addref() { ++refcnt; return this; }
	int release()
	{
		--refcnt;
		if(refcnt == 0) { delete this; return 0; }
		return refcnt;
	}

	// upload/download to/from the device
	void *syncTo(bool device);
	void *syncToDevice() { return syncTo(true); }
	void *syncToHost() { return syncTo(false); }

	// texture access
	void bind_texture(textureReference &texref);
	void unbind_texture(textureReference &texref);

private:
	// garbage collection facilities
	struct allocated_pointers : public std::set<cuxSmartPtr_impl_t *>
	{
		~allocated_pointers();
	};
	static allocated_pointers all_cuxSmartPtrs;
	static void global_gc();

private:
	// constructs a CUDA array. Used by bind_texture
	cudaArray *getCUDAArray(cudaChannelFormatDesc &channelDesc);

	// garbage collection -- release all unused copies of this pointer
	// on devices other than master
	void gc();

	// ensure no shenanigans (no copy constructor & operator)
	cuxSmartPtr_impl_t(const cuxSmartPtr_impl_t &);
	cuxSmartPtr_impl_t& operator=(const cuxSmartPtr_impl_t &);
};

/*
	Smart pointer to an array in GPU or CPU memory, with on-demand garbage collection, and auto-syncing of CPU/GPU copies

		- Points to an arrayPtr<>-styled n-dimensional block of memory
		- Is reference-counted, and automatically releases all memory when
		the reference count falls to zero.
		- At any given time, the block is either on GPU, CPU, or bound to texture (in a CUDA array)
		- To access the block on host/device, obtain a hptr<> or gptr<> via cast.
		Obtaining will moves the block to the host/device if necessary, invalidating
		any previously obtained pointers to device/host. The pointer is then
		"locked" to that particular device.
		- For convenience, data on the host can be accessed directly off this object,
		with the() operator. This is equivalent to obtaining a hptr<>.
		- Blocks remain allocated on device and host (for speed) until an OOM
		condition is reached on the device. Garbage collection will then remove
		all unlocked device pointers, and retry the allocation.

		- Can be bound to a texture using bind_texture. Must be unbound
		before deallocation. Direct use is discouraged; use cuxTextureBinder
		where possible.
*/
template<typename T>
struct cuxSmartPtr
{
protected:
	cuxSmartPtr_impl_t *m_impl;

public:
	//
	// iterator-ish interface (host)
	//
	struct iterator
	{
		cuxSmartPtr<T> *d;
		uint32_t i, j, k;

		iterator(cuxSmartPtr<T> *d_ = NULL, uint32_t i_ = 0, uint32_t j_ = 0, uint32_t k_ = 0)
		: d(d_), i(i_), j(j_), k(k_)
		{}

		T& operator *()  const { return (*d)(i, j, k); }
		T& operator ->() const { return (*d)(i, j, k); }
		iterator &operator ++()	// prefix
		{
			if(++i >= d->width())
			{
				i = 0;
				if(++j >= d->height())
				{
					j = 0;
					++k;
				}
			}
			return *this;
		}
		bool operator==(const iterator &a) const
		{
			return a.i == i && a.j == j && a.k == k;
		}
		bool operator!=(const iterator &a) const
		{
			return !(a == *this);
		}
	};
	iterator begin() { return iterator(this, 0, 0, 0); }
	iterator end() { return iterator(this, 0, 0, depth()); }
public:
	uint32_t elementSize() const { return m_impl->m_elementSize; }		// size of the stored element (typically, sizeof(T))
	uint32_t width()       const { return m_impl->m_width; }		// logical width of the data array
	uint32_t pitch()       const { return m_impl->m_data.extent[0]; }	// byte width of the data array
	uint32_t height()      const { return m_impl->m_data.extent[1]; }	// logical height of the data array
	uint32_t depth()       const { return m_impl->m_data.extent[2]; }	// logical depth of the data array
	uint32_t size()        const { return width()*height()*depth(); }	// logical size (number of elements) in the data array
	uint32_t memsize()     const { return m_impl->memsize(); }		// physical size, in bytes, of the array (including any padding)
	uint32_t extent(int n) const						// logical size (number of elements) of the requested dimension (0=first, 1=second, ...)
		{ return n == 0 ? width() : m_impl->m_data.extent[n]; }

	cuxSmartPtr(const cuxSmartPtr& t)				// construct a copy of existing pointer
	{
		m_impl = t.m_impl->addref();
	}
	cuxSmartPtr &operator =(const cuxSmartPtr &t)			// construct a copy of existing pointer
	{
		if(t.m_impl != m_impl)
		{
			m_impl->release();
			m_impl = t.m_impl->addref();
		}
		return *this;
	}
	~cuxSmartPtr()
	{
		m_impl->release();
	}

	cuxSmartPtr<T> clone(bool copyData = false) const	// make an exact copy of this pointer (the array shape), optionally copying the data as well
	{
		cuxSmartPtr<T> ret(new cuxSmartPtr_impl_t(elementSize(), pitch(), width(), height(), depth()));
		if(copyData)
		{
			syncToHost();
			ret.syncToHost();
			memcpy(ret.m_impl->m_data.ptr, m_impl->m_data.ptr, m_impl->memsize());
		}
		return ret;
	}

	// comparisons/tests
	operator bool() const			// return true if this is not a null pointer
	{
		return (bool)(*m_impl);
	}

	// host data accessors (note: use hptr<> interface if maximum speed is needed)
	T& operator()(const uint32_t x, const uint32_t y, const uint32_t z) const	// 3D accessor (valid only if dim = 3)
	{
		assert(x < width());
		assert(y < height());
		assert(z < depth());
		if(m_impl->onDevice || !m_impl->m_data.ptr) { syncToHost(); }
		return (*(arrayPtr<T, 3> *)(&m_impl->m_data))(x, y, z);
	}
	T& operator()(const uint32_t x, const uint32_t y) const	// 2D accessor (valid only if dim >= 2)
	{
		assert(x < width());
		assert(y < height());
		if(m_impl->onDevice || !m_impl->m_data.ptr) { syncToHost(); }
		return (*(arrayPtr<T, 2> *)(&m_impl->m_data))(x, y);
	}
	T& operator()(const uint32_t x) const			// 1D accessor (valid only if dim >= 2)
	{
		assert(x < width());
		if(m_impl->onDevice || !m_impl->m_data.ptr) { syncToHost(); }
		return (*(arrayPtr<T, 1> *)(&m_impl->m_data))(x);
	}

	// texture access
	void bind_texture(textureReference &texref)
	{
		m_impl->bind_texture(texref);
	}
	void unbind_texture(textureReference &texref)
	{
		m_impl->unbind_texture(texref);
	}
public:
	template<int dim>
		operator hptr<T, dim>()			// request access to data on the host
		{
			syncToHost();
			return *(hptr<T, dim> *)(&m_impl->m_data);
		}

	template<int dim>
		operator gptr<T, dim>()			// request access to data on the GPU
		{
			syncToDevice();
			return *(gptr<T, dim> *)(&m_impl->m_data);
		}

	cuxSmartPtr(uint32_t width = 0, uint32_t height = 1, uint32_t depth = 1, uint32_t elemSize = 0xffffffff, uint32_t align = 128)
	{
		if(elemSize == 0xffffffff) { elemSize = sizeof(T); }
		m_impl = new cuxSmartPtr_impl_t(elemSize, roundUpModulo(elemSize*width, align), width, height, depth);
	}

	void reset()	// clear the value of this pointer
	{
		*this = cuxSmartPtr();
	}

protected:
	// multi-device support. Don't call these directly; use gptr instead.
	T* syncToHost() const { return (T*)const_cast<cuxSmartPtr<T> *>(this)->m_impl->syncToHost(); }
	T* syncToDevice() { return (T*)m_impl->syncToDevice(); }

	// protected constructors
	cuxSmartPtr(cuxSmartPtr_impl_t *d)
	{
		m_impl = d;
	}
};

// copy data from a host C pointer to an cuxSmartPtr<>. Assumes the host array has the
// required number of elements and no stride.
template<typename T>
void copy(cuxSmartPtr<T> &p, const T* data)
{
	if(data == NULL) { return; }

	hptr<T, 3> v = p;
	for(int i = 0; i != p.width(); i++)
	{
		for(int j = 0; j != p.height(); j++)
		{
			for(int k = 0; k != p.depth(); k++)
			{
				v(i, j, k) = *data;
				data++;
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////
// Texturing support
////////////////////////////////////////////////////////////////////////////////////////////////

/*
	Work around CUDA defficiency with some built-in struct alignments.

	CUDA header files declare some structs (float2 being an example) with
	__builtin_align() attribute that resolves to __align__ only when using
	CUDACC. This makes those structure's alignments different in nvcc compiled
	object code, compared to GCC's. Example: with nvcc, float2 is 8-byte
	aligned; on gcc, it's 4-byte aligned (given its members are all 4-byte
	aligned floats). Therefore, a struct that has a float2 member may be
	packed differently on gcc and nvcc. Example: struct { float a; float2 b; };
	On nvcc, &b = &a + 8 (in bytes). On gcc, &b = &a + 4 (bytes).

	This cludge works around the problem by deriving an aligned type from
	the problematic CUDA type. It should be used instead of the CUDA type
	in structures where this problem may occur.
*/
struct ALIGN(8) afloat2 : public float2
{
	afloat2& operator =(const float2 &a) { (float2&)*this = a; return *this; }
};

/*
	cuxTexture<T,dim> -- host copy of texture data and associated texture coordinates

	Convenience class -- holds both the pointer to texture data
	and the texture coordinates in a same object.
*/
template<typename T, int dim=1>
struct cuxTexture : public cuxSmartPtr<T>
{
	// texture coordinates
	afloat2		coords[dim];

	// default and copy constructor
	cuxTexture()						{ }
	cuxTexture(const cuxTexture<T, dim>& a)			{ *this = a; }

	// other commonly used constructors
	explicit cuxTexture(const cuxSmartPtr<T>& a)					 { set(a, make_float2(0., 1.)); }
	explicit cuxTexture(const cuxSmartPtr<T>& a, float2 tc)				 { set(a, tc); }
	explicit cuxTexture(const cuxSmartPtr<T>& a, const float2 *tc)			 { set(a, tc); }
	explicit cuxTexture(const cuxSmartPtr<T>& a, float2 tcx, float2 tcy)		 { set(a, tcx, tcy); }
	explicit cuxTexture(const cuxSmartPtr<T>& a, float2 tcx, float2 tcy, float2 tcz) { set(a, tcx, tcy, tcz); }

	explicit cuxTexture(size_t nx, float2 tcx = make_float2(0., 1.))
		{ set(cuxSmartPtr<T>(nx), tcx); }
	explicit cuxTexture(size_t nx, size_t ny, float2 tcx = make_float2(0., 1.), float2 tcy = make_float2(0., 1.))
		{ set(cuxSmartPtr<T>(nx, ny), tcx, tcy); }
	explicit cuxTexture(size_t nx, size_t ny, size_t nz, float2 tcx = make_float2(0., 1.), float2 tcy = make_float2(0., 1.), float2 tcz = make_float2(0., 1.))
		{ set(cuxSmartPtr<T>(nx, ny, nz), tcx, tcy, tcz); }

	// assignment
	cuxTexture<T, dim>& operator =(const cuxTexture<T, dim>& a)	{ set(a, a.coords); return *this; }
	cuxTexture<T, dim>& operator =(const cuxSmartPtr<T>& a)		{ (cuxSmartPtr<T> &)(*this) = a; return *this; }

	// texture initialization/setting
	void set(const cuxSmartPtr<T>& a, const float2 tc)	// NOTE: This method initializes _all_ texture coordinates to tc
	{
		*this = a;
		for(int i=0; i != dim; i++) { coords[i] = tc; }
	}
	void set(const cuxSmartPtr<T>& a, const float2 tcx, const float2 tcy)
	{
		assert(dim == 2);

		*this = a;
		coords[0] = tcx;
		coords[1] = tcy;
	}
	void set(const cuxSmartPtr<T>& a, const float2 tcx, const float2 tcy, const float2 tcz)
	{
		assert(dim == 3);

		*this = a;
		coords[0] = tcx;
		coords[1] = tcy;
		coords[2] = tcz;
	}
	void set(const cuxSmartPtr<T>& a, const float2 *tc)
	{
		*this = a;
		for(int i=0; i != dim; i++) { coords[i] = tc[i]; }
	}
};

/*
	cuxTextureReferenceInterface -- abstract bind/unbind interface for textures

	Common, type-independent base class for textures, permitting cuxTextureBinder
	to bind/unbind the texture without the need to know the texture type.

	Base class of cuxTextureReference.
*/
struct cuxTextureReferenceInterface
{
	virtual void bind(const void *data_, const float2 *texcoord) = 0;
	virtual void unbind() = 0;
};

// CPU emulation & management
template<int dim>
	struct cuxTexCoords
	{
		afloat2 tc[dim];
		__host__ const float2 &operator[](int i) const { return tc[i]; }
		__host__       float2 &operator[](int i)       { return tc[i]; }
	};

/*
	cuxTextureReference<T,dim,readmode> -- host texture reference interface

	Holds the pointer to CUDA texture reference, and implements host emulation
	of texture sampling. There is a 1-to-1 relationship between a cuxTextureReference<>
	object and a CUDA texture refrence.

	Use bind() and unbind() methods to bind a cuxSmartPtr/cuxTexture to the
	texture reference (or, preferably, cuxTextureBinder).

	You should NEVER have to explicitly instantiate a cuxTextureReference<> object;
	use DEFINE_TEXTURE and DECLARE_TEXTURE macros instead.
*/
template<typename T, int dim, enum cudaTextureReadMode mode>
	struct cuxTextureReference : public cuxTextureReferenceInterface
	{
	public:
/*		cuxSmartPtr<T>	data;
		cuxTexCoords<dim> tc;*/
		cuxTexture<T, dim> tex;		// the texture from which we're going to sample
		textureReference &texref;	// the CUDA texture reference to which the current texture will be bound
		const char *tcSymbolName;	// the symbol name of the __constant__ variable that holds the texture
						// coordinates on device
	public:
		cuxTextureReference(textureReference &texref_, const char *tcSymbolName)
			: texref(texref_), tcSymbolName(tcSymbolName)
		{
		}

		virtual void bind(const void *cuxSmartPtr_data_, const float2 *texcoord)
		{
			bind(*(const cuxSmartPtr<T> *)cuxSmartPtr_data_, texcoord);
		}

		void bind(const cuxSmartPtr<T> &img, const float2 *texcoord)
		{
			unbind();

			tex.set(img, texcoord);

			cuxUploadConst(tcSymbolName, tex.coords);
			tex.bind_texture(texref);
		}

		void bind(const cuxTexture<T, dim> &tex)
		{
			bind(tex, tex.coords);
		}

		virtual void unbind()
		{
			if(!tex) { return; }

			tex.unbind_texture(texref);
			tex.reset();
		}

		float clamp(float i, int d) const
		{
			if(i < 0.f) { i = 0.f; }
			else if(i >= tex.extent(d)) { i = tex.extent(d)-1; }

			return i;
		}

		T tex1D(float x) const		// sample a 1D texture
		{
			// FIXME: implement interpolation, clamp modes, normalized coordinates
			uint32_t i = (uint32_t)clamp(x, 0);

			return tex(i);
		}

		T tex2D(float x, float y) const		// sample from 3D texture
		{
			// FIXME: implement interpolation, clamp modes, normalized coordinates
			uint32_t i = (uint32_t)clamp(x, 0);
			uint32_t j = (uint32_t)clamp(y, 1);

			return tex(i, j);
		}
		T tex3D(float x, float y, float z) const		// sample from 3D texture
		{
			// FIXME: implement interpolation, clamp modes, normalized coordinates
			uint32_t i = (uint32_t)clamp(x, 0);
			uint32_t j = (uint32_t)clamp(y, 1);
			uint32_t k = (uint32_t)clamp(z, 2);

			return tex(i, j, k);
		}
};

/*
	tex1D/2D/3D -- Host implementation of texture sampling, source-compat w. CUDA

	These are here primarely for source-compabitbility with CUDA tex?D calls,
	so that the same kernel code can be recompiled on the host without
	modifications. They simply forward the call to apropriate cuxTextureReference
	methods.
*/
template<typename T, int dim, enum cudaTextureReadMode mode>
	inline T tex1D(const cuxTextureReference<T, dim, mode> &texref, float x)
	{
		return texref.tex1D(x);
	}

template<typename T, int dim, enum cudaTextureReadMode mode>
	inline T tex2D(const cuxTextureReference<T, dim, mode> &texref, float x, float y)
	{
		return texref.tex2D(x, y);
	}

template<typename T, int dim, enum cudaTextureReadMode mode>
	inline T tex3D(const cuxTextureReference<T, dim, mode> &texref, float x, float y, float z)
	{
		return texref.tex3D(x, y, z);
	}

/*
	sample_impl -- GPU implementation of sampling with texture coordinates

	Transforms the "real space" coordinates (x,y,z), to texture-space
	coordinates, based on the offset and scale given in tc, and samples
	the texture.

	Must not be called directly, but through TEX?D() macros (below).
*/
// Sampler routines
template<typename T, typename Texref>
	inline __device__ T sample_impl(Texref &r, float x, float2 tc)
	{
		float xi = (x - tc.x) * tc.y + 0.5f;
		T v = tex1D(r, xi);
		return v;
	}

template<typename T, typename Texref>
	inline __device__ T sample_impl(Texref &r, float x, float y, float2 tcx, float2 tcy)
	{
		float xi = (x - tcx.x) * tcx.y + 0.5f;
		float yi = (y - tcy.x) * tcy.y + 0.5f;
		T v = tex2D(r, xi, yi);
		return v;
	}

template<typename T, typename Texref>
	inline __device__ T sample_impl(Texref &r, float x, float y, float z, float2 tcx, float2 tcy, float2 tcz)
	{
		float xi = (x - tcx.x) * tcx.y + 0.5f;
		float yi = (y - tcy.x) * tcy.y + 0.5f;
		float zi = (z - tcz.x) * tcz.y + 0.5f;
		T v = tex3D(r, xi, yi, zi);
		return v;
	}

/*
	Macros to declare, define and sample from CUDA texture references, on device or the host.

	DEFINE_TEXTURE(name...) -- if compiled with nvcc, defines a CUDA textureReference,
	__constant__ storage for texture coordinates, and a host cuxTextureReference<> object named
	'name##Manager'. The host object handles host-based texture sampling, and
	binding/unbinding of textures to the CUDA texture reference.
	If compiled with gcc, declares an extern reference to the name##Manager object.

	DECLARE_TEXTURE(name...) -- Declares an extern reference to the name##Manager
	object (see DEFINE_TEXTURE).

	TEX?D() -- samples from a textures, taking into account the texture coordinates
	bound with the texture

	Sample usage:
	============
		// in .cu file (defining and sampling)
		DEFINE_TEXTURE(mytexture, float, 3, cudaReadModeElementType, false, cudaFilterModeLinear, cudaAddressModeClamp);

		__global__ somekernel(...)
		{
			float x, y, z;
			...
			float val = TEX3D(mytexture, x, y, z)
		}

		// in .cpp file (binding)
		DECLARE_TEXTURE(mytexture, float, 3, cudaReadModeElementType);

		void somefunction(cuxTexture<float, 3> &tex)
		{
			cuxTextureBinder tb_mytexture(mytexture, tex);
			... call CUDA kernel ...
		}

*/
#if !__CUDACC__
	#define DEFINE_TEXTURE(name, T, dim, mode, norm, fMode, aMode) \
		extern cuxTextureReference<T, dim, mode> name##Manager; \
		static cuxTextureReference<T, dim, mode> &name = name##Manager

	#define DECLARE_TEXTURE(name, T, dim, mode) \
		DEFINE_TEXTURE(name, T, dim, mode, dummy, dummy, dummy)

	#define TEX1D(name, x)       sample(name, x)
	#define TEX2D(name, x, y)    sample(name, x, y)
	#define TEX3D(name, x, y, z) sample(name, x, y, z)

	template<typename T, enum cudaTextureReadMode mode>
		inline __device__ T sample(cuxTextureReference<T, 1, mode> &r, float x)
		{
			return sample_impl<T, cuxTextureReference<T, 1, mode> >(r, x, r.tex.coords[0]);
		}

	template<typename T, enum cudaTextureReadMode mode>
		inline __device__ T sample(cuxTextureReference<T, 2, mode> &r, float x, float y)
		{
			return sample_impl<T, cuxTextureReference<T, 2, mode> >(r, x, y, r.tex.coords[0], r.tex.coords[1]);
		}

	template<typename T, enum cudaTextureReadMode mode>
		inline __device__ T sample(cuxTextureReference<T, 3, mode> &r, float x, float y, float z)
		{
			return sample_impl<T, cuxTextureReference<T, 3, mode> >(r, x, y, z, r.tex.coords[0], r.tex.coords[1], r.tex.coords[2]);
		}
#else
	// real thing
	#define DEFINE_TEXTURE(name, T, dim, mode, norm, fMode, aMode) \
		texture<T, dim, mode> name(norm, fMode, aMode); \
		__constant__ cuxTexCoords<dim> name##TC; \
		cuxTextureReference<T, dim, mode> name##Manager(name, #name "TC")

	#define TEX1D(name, x)       sample(name, x, name##TC)
	#define TEX2D(name, x, y)    sample(name, x, y, name##TC)
	#define TEX3D(name, x, y, z) sample(name, x, y, z, name##TC)

	template<typename T, enum cudaTextureReadMode mode>
		inline __device__ T sample(texture<T, 1, mode> r, float x, float2 tc)
		{
			return sample_impl<T, texture<T, 1, mode> >(r, x, tc);
		}
	template<typename T, enum cudaTextureReadMode mode>
		inline __device__ T sample(texture<T, 2, mode> r, float x, float y, float2 tcx, float2 tcy)
		{
			return sample_impl<T, texture<T, 2, mode> >(r, x, y, tcx, tcy);
		}
	template<typename T, enum cudaTextureReadMode mode>
		inline __device__ T sample(texture<T, 3, mode> r, float x, float y, float z, float2 tcx, float2 tcy, float2 tcz)
		{
			return sample_impl<T, texture<T, 3, mode> >(r, x, y, z, tcx, tcy, tcz);
		}


	template<typename T, enum cudaTextureReadMode mode>
		inline __device__ T sample(texture<T, 1, mode> r, float x, cuxTexCoords<1> tc)
		{
			return sample(r, x, tc[0]);
		}
	template<typename T, enum cudaTextureReadMode mode>
		inline __device__ T sample(texture<T, 2, mode> r, float x, float y, cuxTexCoords<2> tc)
		{
			return sample(r, x, y, tc[0], tc[1]);
		}
	template<typename T, enum cudaTextureReadMode mode>
		inline __device__ T sample(texture<T, 3, mode> r, float x, float y, float z, cuxTexCoords<3> tc)
		{
			return sample(r, x, y, z, tc[0], tc[1], tc[2]);
		}
#endif

/*
	texture loading/construction utilities
*/
cuxTexture<float, 1> construct_texture_by_resampling_1D	(double *X, double *Y, int ndata, int nsamp = 1024);
cuxTexture<float, 1> load_constant_texture_1D			(float val, float x0, float x1);
cuxTexture<float, 3> load_constant_texture_3D(
	float val,
	float x0, float x1,	// range of x texture coordinates [x0,x1]->[0,nx-1]
	float y0, float y1,
	float z0, float z1
);
cuxTexture<float, 1> load_and_resample_texture_1D		(const char *fn, int nsamp = 1024);
cuxTexture<float, 1> load_resample_and_clip_texture_1D		(const char *fn, float clipValue = 0.f, int nsamp = 1024);

float2		     texcoord_from_range			(float imgx, float imgy, float x, float y);

/*
	cuxTextureBinder -- RAII style texture binding/unbinding

	Use the cuxTextureBinder to binds a cuxTexture<> to a cuxTextureReference<>. The object's
	constructor will upload the texture coordinates, and bind the texture. The destructor will
	automatically unbind the texture.

	It is recommended to use this mechanism to bind a texture just before a kernel call.
*/
struct cuxTextureBinder
{
	cuxTextureReferenceInterface &tex;

	template<typename T>
	cuxTextureBinder(cuxTextureReferenceInterface &tex_, const cuxSmartPtr<T> &data, const float2 *texcoord)
		: tex(tex_)
	{
		tex.bind(&data, texcoord);
	}

	template<typename T, int dim>
	cuxTextureBinder(cuxTextureReferenceInterface &tex_, const cuxTexture<T, dim> &tptr)
		: tex(tex_)
	{
		tex.bind(&tptr, tptr.coords);
	}

	~cuxTextureBinder()
	{
		tex.unbind();
	}
};

#endif // cuda_cux__
