/***************************************************************************
 *   Copyright (C) 2005 by Mario Juric   *
 *   mjuric@astro.Princeton.EDU   *
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
#ifndef __astro_system_memorymap_h
#define __astro_system_memorymap_h

#ifdef __linux__
#include <sys/mman.h>
#else
#warning "Using fake memory mapping"
#include "fakemmap.h"
#define mmap fakemmap
#define munmap fakemunmap
#define msync fakemsync
#endif

#include <fcntl.h>
#include <string>
#include <iosfwd>

#include <map>
#include <deque>
#include <list>
#include <vector>
#include <iostream>
#include <stdexcept>

namespace peyton {

namespace system {

class MemoryMap
{
protected:
	std::string filename;
	int fd;

	size_t length;
	void *map;

	bool closefd;
public:
	enum {
		ro = PROT_READ,
		wo = PROT_WRITE,
		rw = PROT_READ | PROT_WRITE,
		none = PROT_NONE,
		exec = PROT_EXEC
	};
	
	enum {
		shared = MAP_SHARED,
		priv = MAP_PRIVATE
	};
	
	static const int pagesize;
	static void pagesizealign(std::ostream &out);
	static void pagesizealign(std::istream &in);
public:
	void open(int fd, size_t length_, size_t offset, int mode, int mapstyle, bool closefd = false);
public:
	MemoryMap();
	MemoryMap(const std::string &filename, size_t length = 0, size_t offset = 0, int mode = ro, int map = shared);

	void open(const std::string &filename, size_t length = 0, size_t offset = 0, int mode = ro, int map = shared);
	void sync();
	void close();

	~MemoryMap();

	operator void *() { return map; }

	size_t size() const { return length; }
};

template<typename T>
class MemoryMapVector : public MemoryMap
{
public:
	typedef T* iterator;
	size_t siz;
public:
	MemoryMapVector() : MemoryMap(), siz(0) {}
	void open(const std::string &filename, size_t size = -1, size_t offset = 0, int mode = ro, int mapstyle = shared)
	{
		MemoryMap::open(filename, sizeof(T)*size, offset, mode, mapstyle);
		if(size < 0) { siz = length/sizeof(T); } else { siz = size; }
	}

	const T &operator[](int i) const { return ((const T *)map)[i]; }
	T &operator[](int i) { return ((T *)map)[i]; }

	iterator begin() { return (T *)map; }
	iterator end() { return ((T *)map) + siz; }

	T& front() { return *(T *)map; }
	T& back() { return ((T *)map) + (siz-1); }

	size_t size() const { return siz; }
	size_t allocated() { return length / sizeof(T); }
};

struct MemoryMapError : public std::runtime_error {
	MemoryMapError(const std::string &msg): std::runtime_error(msg) {};
	virtual ~MemoryMapError() throw() {};
};


} // namespace system
} // namespace peyton

using namespace peyton::system;

#define __peyton_system peyton::system
#endif

