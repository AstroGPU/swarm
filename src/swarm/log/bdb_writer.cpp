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

/*! \file bdb_writer.cpp
 *    \brief Writer plugin to output the log into a Berkeley-DB database.
 *
 *
 */

#include "../common.hpp"

#include "../types/config.hpp"
#include "../plugin.hpp"

#include "log.hpp"
#include "writer.h"

#include "bdb_database.hpp"

namespace swarm { namespace log {

using namespace gpulog;

/**
 *  \brief Writer plugin to output the log into a Berkeley-DB database.
 *
 *  To use this, add following lines to your integration configuration file
 *
 *  log_output = bdb
 *  log_output_db = <fileName>
 *
 *  Replace <fileName> with the name of the output file without extension. the db extension will be added
 *  automatically.
 *
 *
 *
 *  *EXPERIMENTAL*: This class is not thoroughly tested.
 *  \ingroup experimental
 * 
 *
 */
class bdb_writer : public writer
{
	// To do this, first we have to include the BDB in the CMake files
	// in process, we should make a gpulog::ilogstream out of the data
	// and read all the gpulog records and put them in the database
	// now the problem is what are we going to use for the key. Maybe we
	// can use a very simple key and later on we can do secondary databases
	// to make the indexes based on time and system id.

    idx_t current_recno;
    bdb_database db;
//! constructor for bdb_writer
public:
	bdb_writer(const config& cfg):current_recno(1){
		std::string fileName = cfg.require("log_output_db",std::string());
		db.createEmpty(fileName);
	}

        //! Process the log data and put them in the database
	virtual void process(const char *log_data, size_t length) {
		ilogstream stream(log_data,length);
		while(logrecord lr = stream.next()){
            db.put(lr, current_recno);
            current_recno += 1;
		}
	}

        //! Destructor
	~bdb_writer(){
		db.close();
	}
};

//! Initialize the database writer plugin
writer_plugin_initializer< bdb_writer >
	bdb_writer_plugin("bdb", "This is the Berkeley DB writer");

} }
