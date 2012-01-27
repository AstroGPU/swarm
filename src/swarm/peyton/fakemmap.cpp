#include "fakemmap.h"
#include <map>
#include <assert.h>
#include <algorithm>
#include <stdio.h>

using std::min;

struct mapping_t {
	int fd;
	off_t offset;
};

#undef SSIZE_MAX
const size_t SSIZE_MAX = 10240UL;

std::map<size_t,mapping_t> mappings;

void *fakemmap(void *addr, size_t length, int prot, int flags, int fd, off_t offset){
	// Allocate memory for the mapping
	if(addr==0) addr = malloc(length);

	// Set up the file and memory
	char* caddr = (char*)addr;
	lseek(fd, offset, SEEK_SET);

	// Read the contents of the file into memory
	int l = 0;
	for(int i = 0; i < 10000; i++){
	  int s = read(fd,caddr+l,min(length-l,SSIZE_MAX)); 
	  //int s = read(fd,caddr+l,min(length-l,10240UL)); 
		if(s==-1){
			return MAP_FAILED;
		}else if(s==0)
			break;
		else
			l+=s;
	}

	// Save mapping information
	mappings[(size_t)addr].fd = fd;
	mappings[(size_t)addr].offset = offset;

	// Same behavior as mmap
	return addr;
}

int fakemsync(void *addr, size_t length, int flags) {
	// Load mapping information
	int fd = mappings[(size_t)addr].fd;
	off_t offset  = mappings[(size_t)addr].offset;

	// Set up file and memory
	char* caddr = (char*)addr;
	lseek(fd, offset, SEEK_SET);

	// Write memory to file
	int l = 0;
	for(int i = 0; i < 10000; i++){
	  int s = write(fd,caddr+l,min(length-l,SSIZE_MAX)); 
	  //int s = write(fd,caddr+l,min(length-l,10240UL)); 
		if(s==-1)
			return -1;
		else if( s == 0)
			break;
		else
			l+=s;
	}

	// Same behavior
	return 0;
}

int fakemunmap(void *addr, size_t length){
	int r = fakemsync(addr,length,0);
	free(addr);
	return r;
}

