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

#pragma once

#include <cuda_runtime_api.h>


namespace swarm {
namespace hp {


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

template< class T >
struct DeviceAllocator {
	typedef T Elem;
	static void free(T * p) { cudaFree(p); }
	static T *  alloc(size_t s) {  
		T* p;
		cudaMalloc(&p, s * sizeof(T) );
		return p;
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

template< class T >
struct MappedHostAllocator : public HostAllocator<T> {
	static T *  alloc(size_t s) {  
		T* p;
		cudaHostAlloc(&p, s * sizeof(T), cudaHostAllocMapped );
		return p;
	}

	static T* getDevicePointer(T* hp){
		T* dp;
		cudaGetDevicePointer(&dp,hp,0);
		return dp;
	}
};

template< class A, class B>
class AllocPair {
};

template< class A> 
class AllocPair < A, A > {
	typedef typename A::Elem T; 
	static void copy(T* begin, T* end, T* dst){
		A::copy(begin,end,dst);
	}
};

template< class T>
void alloc_copy(DefaultAllocator<T>,DeviceAllocator<T>, T* begin, T* end, T* dst){
	cudaMemcpy(dst, begin, (end-begin)*sizeof(T), cudaMemcpyHostToDevice);
}

template< class T>
void alloc_copy(DeviceAllocator<T>,DefaultAllocator<T>, T* begin, T* end, T* dst){
	cudaMemcpy(dst, begin, (end-begin)*sizeof(T), cudaMemcpyDeviceToHost);
}

/*
template< class Alloc, class OtherAlloc >
void allocator_copy(typename Alloc::Elem* begin, typename Alloc::Elem* end, typename Alloc::Elem* dst){ 
	assert(false); 
}

template< class Alloc>
void allocator_copy< Alloc, Alloc >(typename Alloc::Elem* begin,typename Alloc::Elem* end,typename Alloc::Elem* dst) {
	Alloc::copy(begin,end,dst);
}

template<class T>
void allocator_copy< DefaultAllocator<T>, DeviceAllocator<T> >(T* begin, T* end, T* dst) {
	cudaMemcpy(dst, begin, (end-begin)*sizeof(T), cudaMemcpyHostToDevice);
}

template<class T>
void allocator_copy< DefaultAllocator<T>, HostAllocator<T> >(T* begin, T* end, T* dst) {
	cudaMemcpy(dst, begin, (end-begin)*sizeof(T), cudaMemcpyHostToHost);
}

template<class T>
void allocator_copy< HostAllocator<T>, DeviceAllocator<T> >(T* begin, T* end, T* dst) {
	cudaMemcpy(dst, begin, (end-begin)*sizeof(T), cudaMemcpyHostToDevice);
}

template<class T>
void allocator_copy< DeviceAllocator<T>, HostAllocator<T> >(T* begin, T* end, T* dst) {
	cudaMemcpy(dst, begin, (end-begin)*sizeof(T), cudaMemcpyDeviceToHost);
}
*/

}
}
