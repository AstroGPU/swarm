/*************************************************************************
 * Copyright (C) 2010 by Saleh Dindar and the Swarm-NG Development Team  *
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

/*! \file storage.hpp
 *  \brief Classes to provide selective storage space for arrays on GPU.
 *
*/

/// The main namespace for the Swarm-NG library
namespace swarm {

struct LocalMemory{};
struct SharedMemory{};

template< class T, class Storage >
class storage {
	T* p;
public:
	T& operator* (){
		return *p;
	}
	T* operator-> (){
		return p;
	}
	T& operator[] (const int& i){
		return p[i];
	}
	const static int local_memory_usage = 0;
	const static int shared_memory_usage = 0;
	storage(void* local_memory_pointer, void* shared_memory_pointer);
};

template< class T>
class storage<T, LocalMemory>{
	T* p;
public:
	T& operator* (){
		return *p;
	}
	T* operator-> (){
		return p;
	}
	T& operator[] (const int& i){
		return p[i];
	}
	const static int local_memory_usage = sizeof(T);
	const static int shared_memory_usage = 0;
	storage(void* local_memory_pointer, void* shared_memory_pointer)
		:p((T*) local_memory_pointer){}
	
};

template< class T>
class storage<T, SharedMemory>{
	T* p;
public:
	T& operator* (){
		return *p;
	}
	T* operator-> (){
		return p;
	}
	T& operator[] (const int& i){
		return p[i];
	}
	const static int local_memory_usage =  0;
	const static int shared_memory_usage = sizeof(T);
	storage(void* local_memory_pointer, void* shared_memory_pointer)
		:p((T*) shared_memory_pointer){}
	
};

}
