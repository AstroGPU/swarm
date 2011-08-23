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
 * It allows multiple threads access consecutive elements of an structure
 * by grouping the structure elements together. However, it requires that all
 * the elements of the structure should be of the same size.
 *
 * The Item should provide WARPSIZE and scalar_t
 * for offset calculation
 */
template<class Item, typename _Scalar = typename Item::scalar_t, int _WARPSIZE = Item::WARPSIZE >
struct CoalescedStructArray {
public:
	static const int WARPSIZE = _WARPSIZE;
	typedef Item* PItem;
	typedef _Scalar scalar_t;

private:
	PItem _array;
	size_t _block_count;

public:
	GENERIC CoalescedStructArray () {};
	GENERIC CoalescedStructArray(PItem array, size_t block_count)
		:_array(array),_block_count(block_count){}
	GENERIC Item& operator[] ( const int & i ) {
		size_t block_idx = i / WARPSIZE;
		size_t idx = i % WARPSIZE;
		scalar_t * blockaddr = (scalar_t*) (get() + block_idx);
		return * (Item *) ( blockaddr + idx );
	}
	GENERIC int block_count()const{
		return _block_count ;
	}
	GENERIC int size()const{
		return _block_count * WARPSIZE;
	}
	GENERIC Item * get() {
		return _array;
	}

	GENERIC Item * begin() {
		return get();
	}

	GENERIC Item * end() {
		return get() + _block_count;
	}
};


}
