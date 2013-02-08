/*
 * binary_reader.hpp
 *
 *  Created on: Dec 14, 2012
 *      Author: saleh
 */
#pragma once

#include "swarm/common.hpp"
#include "swarm/log/io.hpp"
#include <istream>

class binary_reader {

	std::istream& _input;

    char* current;
    char* buffer_begin;
    char* buffer_end;

public:
	binary_reader(std::istream& input);



	/*!
	 * Check that header of the file matches what we expect
	 * should be called before we start
	 */
	bool validate();


	/*!
	 * Retrieve the next logrecord, may require loading the next chunk of the file.
     *
     * If there is the end of file, we return an invalid logrecord.
	 */
	gpulog::logrecord next();

    bool ensure(const size_t& len);
    char* readBytes(const size_t& len);

    void readNextChunk();
    ptrdiff_t tellg();
    void seek(ptrdiff_t absolute_position);
    void readChunk(ptrdiff_t current_offset);
};
