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

/*! \file binary_writer.cpp
 *    \brief Defines a writer that writes to Berkeley DB in binary form. 
 *
 *
 */

#include "../common.hpp"

#include "../types/config.hpp"
#include "../plugin.hpp"

#include "io.hpp"
#include "writer.h"

namespace swarm {

/**
 * binary_writer plugin for use in
 * io.cpp
 *
 */
class binary_writer : public swarm::writer
{
protected:
	std::auto_ptr<std::ostream> output;
	std::string rawfn, binfn;

public:
	binary_writer(const config &cfg)
	{
		binfn = cfg.at("log_output");
		if(binfn=="")
			ERROR("Expected filename for writer.")
				rawfn = binfn + ".raw";

		output.reset(new std::ofstream(rawfn.c_str()));
		if(!*output)
			ERROR("Could not open '" + rawfn + "' for writing");

		// write header
		swarm::swarm_header fh(UNSORTED_HEADER_FULL);
		output->write((char*)&fh, sizeof(fh));
	}

	/* 
	 * We postpone sorting the output here. We leave the
	 * file unsorted as it is. No one knows who 
	 * will benefit from sorting the file and it is
	 * very inefficient
	 */
	~binary_writer()
	{
		output.reset(NULL);
	}

	virtual void process(const char *log_data, size_t length)
	{
		// TODO: filter out the printfs
		output->write(log_data, length);
	}
};


writer_plugin_initializer< binary_writer >
	binary_writer_plugin("binary", "This is the binary writer");

}

