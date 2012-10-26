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
/*! \file allocators.hpp
 *  \brief Defines routines to allocate memory using different APIs and copy 
 *  between them
 *  
 *  The classes provide a generalized interface for different memory 
 *  allocators. This is similar to what C++ library uses. The difference
 *  is, device memory allocators allocate a memory that is not accessible
 *  to the CPU. So in addition to allocations we need copy routines between
 *  allocators. 
 *
 *  Every allocator provides three basic actions:
 *   - create
 *   - copy
 *   - free
 *
 *  The copy routine inside the allocator only copies inside the same memory
 *  hierarchy. For memory transfers between domains, we need to use the
 *  function template alloc_copy. alloc_copy is specialized for different
 *  combinations of source and destination to use specific CUDA calls for
 *  copying.
 *
 *  Current allocators are:
 *   - C++ new/delete
 *   - CUDA [GPU] device memory 
 *   - CUDA host memory
 *   - CUDA device-mapped host memory
 *
 */

#pragma once

//! Default allocator that uses C++ new/delete
//! This class uses standard C++ routines for allocation
//! and memory manipulation: new[], delete[] and std::copy.
template< class T >
struct DefaultAllocator {
	typedef T Elem;
	static void free(T * p) { delete[] p; }
	static T *  alloc(size_t s) {  return new T[s]; } 

	static void copy( T* begin, T* end, T* dst ) {
		std::copy ( begin, end, dst );
	}

	static T* clone (T* begin, T* end)  { 
		T* p = alloc( end - begin );
		copy( begin, end, p );
		return p;
	}
};

//! CUDA device memory allocator that uses cudaMalloc,cudaMemcpy,cudaFree
//! It creates a pointer that is allocated on the device. The pointer
//! cannot be used by the caller and should only be passed to a CUDA 
//! kernel. The copy uses cudaMemcpy to transfer data between 2 device
//! arrays.
template< class T >
struct DeviceAllocator {
	typedef T Elem;
	static void free(T * p) { cudaFree(p); }
	static T *  alloc(size_t s) {  
		void* p;
		cudaMalloc(&p, s * sizeof(T) );
		return (T*)p;
	}

	static void copy( T* begin, T* end, T* dst ) {
		cudaMemcpy( dst, begin, (end-begin)*sizeof(T), cudaMemcpyDeviceToDevice );
	}

	static T* clone (T* begin, T* end)  { 
		T* p = alloc( end - begin );
		copy( begin, end, p );
		return p;
	}
};


//! CUDA host memory allocator uses cudaMallocHost,cudaMemcpy,cudaFreeHost
//! Host memory allocator is similar to malloc. The pointers point to 
//! memory that can be used by C++. However, CUDA documentation claims that
//! copying to device memory from a CUDA allocated host array is faster than
//! memory allocated using malloc.
template< class T >
struct HostAllocator {
	typedef T Elem;
	static void free(T * p) { cudaFreeHost(p); }
	static T *  alloc(size_t s) {  
		T* p;
		cudaMallocHost(&p, s * sizeof(T) );
		return p;
	}

	static void copy( T* begin, T* end, T* dst ) {
		cudaMemcpy( dst, begin, (end-begin)*sizeof(T), cudaMemcpyHostToHost );
	}

	static T* clone (T* begin, T* end)  { 
		T* p = alloc( end - begin );
		copy( begin, end, p );
		return p;
	}
};

//! CUDA host memory allocator similar to HostAllocator using device mapped memory
//! A Mapped memory is accessible on host and device. However, the pointers are
//! different and this complicated everything. According to CUDA manual, version
//! 4.0 of CUDA SDK uses unified pointers so there in to need to 
//! map the pointer. In that case, The pointer obtained using this allocator
//! can be passed to a kernel.
template< class T >
struct MappedHostAllocator : public HostAllocator<T> {
	static T *  alloc(size_t s) {  
		T* p;
		cudaHostAlloc(&p, s * sizeof(T), cudaHostAllocMapped );
		return p;
	}

	/*!
	 * This is a non-conforming method for getting the Device Pointer
	 * This method is only needed when dealing with device with compute
	 * capabality less than 2.0
	 */
	static T* getDevicePointer(T* hp){
		T* dp;
		cudaGetDevicePointer(&dp,hp,0);
		return dp;
	}
};


//! Simple copy between the same allocator. Uses the copy
//! method of the allocator.
template< class A, class T> 
void alloc_copy(A,A, T* begin, T* end, T* dst){
	A::copy(begin,end,dst);
}

//! Copy from host memory to device memory
template< class T>
void alloc_copy(DefaultAllocator<T>,DeviceAllocator<T>, T* begin, T* end, T* dst){
	cudaMemcpy(dst, begin, (end-begin)*sizeof(T), cudaMemcpyHostToDevice);
}

//! Copy from device memory to host memory
template< class T>
void alloc_copy(DeviceAllocator<T>,DefaultAllocator<T>, T* begin, T* end, T* dst){
	cudaMemcpy(dst, begin, (end-begin)*sizeof(T), cudaMemcpyDeviceToHost);
}

