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

public:
	binary_reader(std::istream& input);



	/*!
	 * Check that header of the file matches what we expect
	 * should be called before we start
	 */
	bool validate();

	/*!
	 * Determine if there is more log records. This function should be side
	 * effect free.
	 */
	bool has_next();

	/*!
	 * Retrieve the next logrecord, may require loading the next chunk of the file.
	 */
	gpulog::logrecord next();

};
