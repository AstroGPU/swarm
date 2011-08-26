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

#ifndef bits_gpulog_log_h__
#define bits_gpulog_log_h__

namespace gpulog
{

	namespace internal
	{

#ifdef __CUDACC__
		__device__ static inline int global_atomicAdd(int *x, int add) {
			return ::atomicAdd(x,add);
		}
#else
		__host__ static inline int global_atomicAdd(int *x, int add) {
			assert(0); // this must not be called from host code.
			return 0;
		}
#endif

	// device internals encapsulation for log_base<> template
	struct dev_internals
	{
		template<typename T>
			__host__ static void alloc(T* &ret, int num = 1)
			{
				cudaMalloc((void **)&ret, num*sizeof(T));
			}

		template<typename T>
			__host__ static const T get(T *ptr)
			{
				T ret;
				cudaMemcpy(&ret, ptr, sizeof(ret), cudaMemcpyDeviceToHost);
				return ret;
			}

		template<typename T>
			__host__ static void set(T *ptr, const T& val)
			{
				cudaMemcpy(ptr, &val, sizeof(*ptr), cudaMemcpyHostToDevice);
			}

		template<typename T>
			__host__ static void dealloc(T* p, int num = 1)
			{
				cudaFree(p);
			}

		__device__ static inline int threadId()
		{
		#ifdef __CUDACC__
			return ((blockIdx.z * gridDim.y + blockIdx.y) * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;
		#else
			assert(0); // this must not be called from host code.
		#endif
		}

#ifdef __CUDACC__
        __device__ static inline int atomicAdd(int *x, int add) {
			return global_atomicAdd(x, add);
		}
#else
        __host__ static inline int atomicAdd(int *x, int add) {
			return global_atomicAdd(x, add);
		}
#endif
	};

	// host internals encapsulation for log_base<> template
	struct host_internals
	{
		template<typename T>
			static void alloc(T* &ret, int num = 1)
			{
				ret = num == 1 ? new T : new T[num];
			}

		template<typename T>
			static const T get(T *ptr)
			{
				return *ptr;
			}

		template<typename T>
			static void set(T *ptr, const T& val)
			{
				*ptr = val;
			}

		template<typename T>
			static void dealloc(T* p, int num = 1)
			{
				if(num == 1) delete p;
				else delete [] p;
			}

		static inline int atomicAdd(int *x, int add) { int tmp = *x; *x += add; return tmp; }
		static int threadId() { return -1; }
	};

	/*
	   workaround for CUDA 2.2 template parsing bug -- CUDA 2.2 tries to compile a template
	   function as __host__ if it returns a derived type T* (or T&, etc..)
	*/
	template<typename T> struct ptr_t
	{	
		T* ptr;
		__host__ __device__ inline ptr_t(T*p) : ptr(p) {}
		__host__ __device__ operator T*() const { return ptr; }
	};
	/* CUDA 2.2 compatible version */
	#define PTR_T(T)  gpulog::internal::ptr_t<T>
	/* CUDA 2.3 and beyond */
	// #define PTR_T(T)  T*

	/*
	   Log template with the write() implementations.
	*/
	template<typename A>
	struct log_base
	{
	protected:
		char *buffer;
		int *at;
		int buf_len;

	public: /* manipulation from host */
		__host__ void alloc(size_t len)
		{
			buf_len = len;

			A::alloc(at, 1);
			A::set(at, 0);
			A::alloc(buffer, len);

			DHOST( std::cerr << "Allocated " << len << " bytes.\n"; )
		}

		__host__ void free()
		{
			A::dealloc(buffer, buf_len);
			A::dealloc(at);
		
			buffer = NULL; buf_len = 0;
			at = NULL;
		}

		__host__ void clear()	/* clear the output buffer (from host side) */
		{
			A::set(at, 0);
		}

		__host__ int fetch_size() const  /* get the size (in bytes) of the data in the output buffer */
		{
			return A::get(at);
		}

		__host__ void set_size(int pos) const  /* move the buffer pointer to position pos */
		{
			A::set(at, pos);
		}

	public: /* manipulation from device */
		__host__ __device__ int capacity() const
		{
			return buf_len;
		}
		
		__host__ __device__ int size() const  /* get the size (in bytes) of the data in the output buffer */
		{
			return *at;
		}

		__host__ __device__ void seek(int pos) const  /* move the buffer pointer to position pos */
		{
			*at = pos;
		}

		__host__ __device__ char* internal_buffer() const  /* return the internal output buffer (mainly for use by copy()) */
		{
			return buffer;
		}

		// test if the buffer has overflowed
		__host__ __device__ inline bool has_overflowed(int idx)
		{
			return idx > buf_len;
		}

	#if 0
		template<typename T1, typename T2, typename T3>
		__device__ inline PTR_T(SCALAR(T3)) write(const int msgid, const T1 &v1, const T2 &v2, const T3 &v3)
		{
			typedef internal::pktsize<header, T1, T2, T3> P;
			P::dump();

			// allocate and test for end-of-buffer
			int len = P::len_with_padding(v3);
			int at = A::atomicAdd(this->at, len);
			if(has_overflowed(at + len)) { return NULL; }
			char *ptr = buffer + at;

			// write
			header v0(msgid, len);
			P::IO0::put(ptr, v0, P::begin0, P::len0);
			P::IO1::put(ptr, v1, P::begin1, P::len1);
			P::IO2::put(ptr, v2, P::begin2, P::len2);
			P::IO3::put(ptr, v3, P::begin3, P::len3);

			#if ARGINFO
			P::store_arginfo(ptr, v3);
			#endif

			DHOST( std::cerr << "Total packet len = " << len << "\n"; )
			return (SCALAR(T3)*)(ptr + P::begin3);
		}
	#else
		#include "gpulog_write.h"
	#endif


	};

	// Device specialization
	typedef log_base<dev_internals> device_log;

	// Host specialization, with memory management
	struct host_log : public log_base<host_internals>
	{
		host_log(size_t len = 0)
		{
			alloc(len);
		}

		~host_log()
		{
			free();
		}
	};


	/*
		Log memory management and copying API
	*/

	// Download the pointers from symbol 'name' to device_log structure
	// on the host
	inline void download_device_log(device_log &log, device_log* dlog)
	{
		cudaMemcpy(&log, dlog, sizeof(log), cudaMemcpyDeviceToHost);
	}

	// Download the pointers from symbol 'name' to device_log structure
	// on the host
	inline void download_device_log(device_log &log, const char *name)
	{
		cudaMemcpyFromSymbol(&log, name, sizeof(log), 0, cudaMemcpyDeviceToHost);
	}

	// Uplaod the pointers from device_log structure on the host to
	// symbol 'name'
	inline void upload_device_log(const char *name, device_log &log)
	{
		cudaMemcpyToSymbol(name, &log, sizeof(log), 0, cudaMemcpyHostToDevice);
	}

	// Uplaod the pointers from device_log structure on the host to a new buffer
	inline device_log* upload_device_log(device_log &log)
	{
		void* pdlog;
		cudaMalloc(&pdlog, sizeof(log));
		cudaMemcpy(pdlog, &log, sizeof(log), cudaMemcpyHostToDevice);
		return (device_log*) pdlog;
	}

	//
	// Copy device log buffers and metadata from 'from'
	// to host_log 'to'.
	//
	// If (flags & LOG_DEVCLEAR) reset the device log afterwards
	//
	__host__ inline void copy(host_log &to, device_log &from, int flags = 0)
	{
		// clear host log
		to.clear();

		// memcpy from device log
		int size = from.fetch_size();
		if(size == 0) { return; }

		// clear/resize host log if needed
		if(to.capacity() != from.capacity())
		{
			to.free();
			to.alloc(from.capacity());
		}

		// memcpy the data
		cudaMemcpy(to.internal_buffer(), from.internal_buffer(), size, cudaMemcpyDeviceToHost);
		to.set_size(size);

		// clear device log if asked for
		if(flags & LOG_DEVCLEAR)
		{
			from.clear();
		}
	}

	//
	// Copy device log buffers and metadata from object stored
	// in device symbol 'symbol', to host_log to.
	//
	// If (flags & LOG_DEVCLEAR) reset the device log afterwards
	//
	inline void copy(host_log &to, const char *from, int flags = 0)
	{
		device_log dlog;
		download_device_log(dlog, from);
		copy(to, dlog, flags);
	}

	//
	// Copy device log buffers and metadata from object stored
	// in device symbol 'symbol', to host_log to.
	//
	// If (flags & LOG_DEVCLEAR) reset the device log afterwards
	//
	inline void copy(host_log &to, device_log *from, int flags = 0)
	{
		device_log dlog;
		download_device_log(dlog, from);
		copy(to, dlog, flags);
	}

	//
	// Setup a new log object in device symbol 'symbol'
	//
	inline device_log alloc_device_log(const char *symbol, size_t len)
	{
		device_log dlog;
		dlog.alloc(len);
		upload_device_log(symbol, dlog);
		return dlog;
	}

	//
	// Setup a new log object and return the device pointer
	//
	inline device_log* alloc_device_log(size_t len)
	{
		device_log dlog;
		dlog.alloc(len);
		return upload_device_log(dlog);
	}

	//
	// Free the memory associated with device log in symbol 'symbol'
	//
	inline void free_device_log(const char *symbol)
	{
		device_log dlog;
		download_device_log(dlog, symbol);
		dlog.free();
	}



	} // namespace internal
} // namespace gpulog

#endif //  bits_gpulog_log_h__
