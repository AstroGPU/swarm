#include "binary_reader.hpp"

#include <exception>

using namespace swarm;
using namespace std;
using gpulog::logrecord;

const int MAX_RECORD_SIZE = 4096;
const int BUFFER_SIZE = 128*1024*1024;
const int PAGESIZE = 4096;

/**
 * We should read from the file into a buffer
 */

binary_reader::binary_reader(istream& input)
    :_input(input)
{
        buffer_end = current = buffer_begin = new char[BUFFER_SIZE];
}

ptrdiff_t binary_reader::tellg(){
    return _input.tellg() - (buffer_end - current);
}

void binary_reader::seek(ptrdiff_t absolute_position){
    ptrdiff_t page_offset = ((absolute_position + PAGESIZE-1)/PAGESIZE)*PAGESIZE;
    ptrdiff_t mode_offset = absolute_position - page_offset;

    _input.seekg( page_offset, ios_base::beg );

    readChunk( mode_offset );
}

void binary_reader::readChunk(ptrdiff_t current_offset){
    // read the chuck
    _input.read(buffer_begin,BUFFER_SIZE);
    buffer_end = buffer_begin + _input.gcount();

    // set the current pointer
    current = buffer_begin + current_offset;

    // std::cerr << "Read " << (buffer_end-buffer_begin) << 
    //   " bytes and current at " << (current-buffer_begin) << std::endl;
}

void binary_reader::readNextChunk(){
    // lseek back to before of the current
    ptrdiff_t back_offset = buffer_end - current;
    ptrdiff_t page_offset = ((back_offset + PAGESIZE-1)/PAGESIZE)*PAGESIZE;
    if(page_offset > 0) _input.seekg( -page_offset , ios_base::cur);

    readChunk(page_offset - back_offset);
}


bool binary_reader::ensure(const size_t& len){
    if(current + len < buffer_end){
        return true;
    } else if(!_input.eof()){
        readNextChunk();
        return true;
    } else {
        return false;
    }
}

char* binary_reader::readBytes(const size_t& len){
    if(ensure(len)){
        char* ret = current;
        current += len;
        return ret;
    } else {
        return 0;
    }
}

bool binary_reader::validate() {
	swarm_header defh(query::UNSORTED_HEADER_FULL);

    readNextChunk();

	swarm_header* curh = reinterpret_cast<swarm_header*>(readBytes(sizeof(swarm_header)));

	return _input.good() && curh && defh.is_compatible(*curh);
}



logrecord binary_reader::next() {
    const static gpulog::internal::header hend(-1,0);
	const static int LOH = sizeof(int)*2;

    if(!ensure(LOH)) return logrecord((char*)&hend);
	logrecord l(current);

	//std::cerr << "Len: " << l.len() << std::endl;
	if(l.len() > MAX_RECORD_SIZE) {
		std::cerr << "Trying to read very large error, probably wrong data: " << l.len() << std::endl;
		throw std::runtime_error("Record is too large");
	}

    // We have to re-read the logrecord since the current might move
    logrecord ll(readBytes(l.len()));

	return ll;
}

