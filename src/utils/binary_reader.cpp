#include "binary_reader.hpp"

#include <exception>

using namespace swarm;
using namespace std;
using gpulog::logrecord;

const int MAX_RECORD_SIZE = 4096;

/**
 * We should read from the file into a buffer
 */

binary_reader::binary_reader(istream& input):_input(input){}

bool binary_reader::validate() {
	swarm_header defh(query::UNSORTED_HEADER_FULL);
	swarm_header curh(query::UNSORTED_HEADER_FULL);

	_input.read(reinterpret_cast<char*>(&curh),sizeof(curh));

	return _input.good() && defh.is_compatible(curh);
}





bool binary_reader::has_next() {
	return !_input.eof();
}


logrecord binary_reader::next() {
	static char buffer[MAX_RECORD_SIZE];
	const static int LOH = sizeof(int)*2/sizeof(char);

	// The buffer is the logrecord, we will read the logrecord into the
	// buffer with the aliased pointers.
	logrecord l(buffer);

	// First read the header into the buffer, this fill in the
	// length of logrecord
	_input.read(buffer, LOH);

	//std::cerr << "Len: " << l.len() << std::endl;
	if(l.len() > MAX_RECORD_SIZE) {
		std::cerr << "Trying to read very large error, probably wrong data: " << l.len() << std::endl;
		throw std::runtime_error("Record is too large");
	}

	// Then read the rest of the record, since we have alread read the header and there is
	// no going back, we just read the rest.
	_input.read(buffer +LOH, l.len() - LOH);

	return l;
}

