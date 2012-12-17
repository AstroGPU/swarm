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

/*! \file coalescedstructarray.hpp
 *  \brief Defines and implements the array of structures with coalescsed access. 
 */

#pragma once


//! \todo these defines should be in a global header
//! they are placed here because we don't have global header
#ifdef __CUDACC__
#define GENERIC inline __device__ __host__
#define GPUAPI inline __device__
#else
#define GENERIC inline
#define GPUAPI inline 
#endif 


namespace swarm {

/*!
 * Array of structures with coalecsed access.
 *
 * It allows multiple threads access consecutive elements of a structure
 * by grouping the structure elements together. However, it requires that all
 * the elements of the structure should be of the same size.
 *
 * The Item should provide CHUNK_SIZE and scalar_t
 * for offset calculation
 *
 * The Item should be designed with care. Item should be a struct
 * where each of its members an array of scalar_t[CHUNK_SIZE]. 
 * That is the granularity level of the Item. The Item however,
 * can have 2D arrays or other arrays. But the scalar_t[CHUNK_SIZE] 
 * is the smallest item. For example, an array of scalar_t[2*CHUNK_SIZE]
 * and a 2D array of scalar_t[3][CHUNK_SIZE]
 * are allowed. Item can contain other structs with the same granularity.
 * The following example shows this feature:
 * struct A {
 *   double _x[8];
 *   double _y[8];
 *   double _z[16];
 *   double& z(const int& i){
 *   	return z[i*8];
 *   }
 *   double& x(){
 *   	return _x[0];
 *   }
 * };
 *
 * struct Item {
 *   typedef double scalar_t;
 *   static const int CHUNK_SIZE = 8;
 *   A _a[4];
 *   double _f[8];
 *   double _m[16];
 *   double& m(const int& i){
 *   	return m[i*8];
 *   }
 *   double& f(){
 *   	return _f[0];
 *   }
 *   A& a(const int& i){
 *   	return a[i];
 *   }
 * }
 *
 * The Item should also have accessors that return the first
 * Item of the arrays and other items should not be accessible.
 * Look at the example above.
 *
 */
template<class Item, typename _Scalar = typename Item::scalar_t, int _CHUNK_SIZE = Item::CHUNK_SIZE >
struct CoalescedStructArray {
public:
	static const int CHUNK_SIZE = _CHUNK_SIZE;
	//! We may use a shared_ptr for PItem
	typedef Item* PItem;
	typedef _Scalar scalar_t;

private:
	PItem _array;
	size_t _block_count;

public:
	GENERIC CoalescedStructArray () {};
	/*!
	 * CoalescedStructArray lets you allocate the memory for
	 * the array. the second argument is called block_count 
	 * because the real number of items in the array is
	 * block_count*CHUNK_SIZE, since there are CHUNK_SIZE items
	 * in each block.
	 * This constructor can be used as:
	 * \code
	 * CoalescedStructArray&lt;Item&gt; arr(new Item[n], n);
	 * \endcode
	 *
	 * Memory management is left to the user because CoalescedStructArray
	 * is meant to be used by GPU functions and memory management in that
	 * case is very complicated.
	 *
	 */
	GENERIC CoalescedStructArray(PItem array, size_t block_count)
		:_array(array),_block_count(block_count){}

	/**
	 * The core functionality: overlapping Items to allow 
	 * transparent coalesced access. To calculate the address of
	 * item i, i is divided by CHUNK_SIZE, the quotient is the 
	 * block address and remainder is the offset into the item.
	 * For example, if the remainder is 3, it means that we want
	 * to access the 3rd item of all scalar_t[CHUNK_SIZE] arrays.
	 * The trick is to slide the pointers by that amount so the 
	 * item that we want to access becomes the first item so we
	 * can use the default accessors on the item.
	 */
	GENERIC Item& operator[] ( const int & i ) {
		size_t block_idx = i / CHUNK_SIZE;
		size_t idx = i % CHUNK_SIZE;
		scalar_t * blockaddr = (scalar_t*) (get() + block_idx);
		return * (Item *) ( blockaddr + idx );
	}

	/**
	 * Number of blocks in this array
	 *
	 */
	GENERIC int block_count()const{
		return _block_count ;
	}
	/**
	 * Number of items in the array: (block-count*CHUNK_SIZE)
	 */
	GENERIC int size()const{
		return _block_count * CHUNK_SIZE;
	}

	/** 
	 * The raw pointer
	 */
	GENERIC Item * get() {
		return _array;
	}

	/**
	 * Begin pointer, it is a raw pointer. It should 
	 * only be used for copying the whole array
	 */
	GENERIC Item * begin() {
		return get();
	}

	/**
	 * End pointer, it is a raw pointer. It should 
	 * only be used for copying the whole array
	 */
	GENERIC Item * end() {
		return get() + _block_count;
	}
};

//! CoalescedStruct for double data type
template<int W>
struct DoubleCoalescedStruct {
	typedef double scalar_t;
	double _value[W];
	GENERIC double& value(){ return _value[0]; }
};

}
