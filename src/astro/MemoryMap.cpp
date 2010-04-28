//
// Author: Mario Juric <mjuric@cfa.harvard.edu>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//

#include <astro/memorymap.h>
#include <astro/util.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <iostream>
#include "swarm.h"

using namespace std;

const int MemoryMap::pagesize = getpagesize();

void MemoryMap::pagesizealign(ostream &out)
{
	size_t at = out.tellp();
	size_t offs = at % pagesize;
	if(offs) { out.seekp(pagesize - offs, ios::cur); }
}

void MemoryMap::pagesizealign(istream &in)
{
	size_t at = in.tellg();
	size_t offs = at % pagesize;
	if(offs) { in.seekg(pagesize - offs, ios::cur); }
}

void MemoryMap::open(const std::string &fn, size_t length_, size_t offset, int mode, int mapstyle)
{
	if(offset > 0 && (offset % pagesize)) { ERROR("Invalid offset requested for memory mapped area - not a multiple of pagesize (" + str(pagesize) + ")"); }

	int flags = 0;

	     if(mode & rw) { flags |= O_RDWR | O_CREAT; }
	else if(mode & ro) { flags |= O_RDONLY; }
	else if(mode & wo) { flags |= O_WRONLY | O_CREAT; }
	else ERROR("Invalid mode parameter - mode needs to include ro, wo or rw");

	int fd = ::open(fn.c_str(), flags);
	if(fd == -1)
	{
		fd = 0;
		ERROR(string("Error opening file [") + fn + "]");
	}

	if(length_ == 0 && !(flags & O_WRONLY))
	{
		struct stat buf;
		fstat(fd, &buf);
		length_ = buf.st_size;
	}

	close();
	filename = fn;
	open(fd, length_, offset, mode, mapstyle, true);
}

void MemoryMap::open(int fd_, size_t length_, size_t offset, int prot, int mapstyle, bool closefd_)
{
	close();
	length = length_;
	closefd = closefd_;
	fd = fd_;

	// check if the length of the file is sufficient to hold the requested length
	// if not, enlarge if possible
	struct stat buf;
	fstat(fd, &buf);
	if(buf.st_size < offset + length)
	{
		if(!(prot & PROT_WRITE))
		{
			close();
			ERROR(string("A read-only file is smaller than the length requested to be memory mapped"));
		}

		// resize the file
		lseek(fd, offset + length - 1, SEEK_SET);
		char b = 1;
		write(fd, &b, 1);
	}

	map = mmap(0, length, prot, mapstyle, fd, offset);
	if(map == MAP_FAILED)
	{
		map = NULL;
		close();
		ERROR(string("Memory mapping of file [") + filename + "] falied. Parameters: length=" + str(length) + ", offset=" + str(offset));
	}
}

void MemoryMap::sync()
{
	assert(fd != 0);
	assert(map != NULL);
	msync(map, length, MS_SYNC);
}

MemoryMap::MemoryMap()
: filename(""), fd(0), map(NULL), length(0), closefd(true)
{
}

MemoryMap::MemoryMap(const std::string &fn, size_t length_, size_t offset, int mode, int mapstyle)
: filename(fn), fd(0), map(NULL), length(length_)
{
	open(fn, length, offset, mode, mapstyle);
}

void MemoryMap::close()
{
	if(map != NULL)
	{
		if(munmap(map, length) == -1) { ERROR(string("Error unmapping file [") + filename + "]"); }
		map = NULL;
	}

	if(closefd && fd != 0)
	{
		if(::close(fd) == -1) { ERROR(string("Error closing file [") + filename + "]"); }
		fd = 0;
	}
}

MemoryMap::~MemoryMap()
{
	close();
}
