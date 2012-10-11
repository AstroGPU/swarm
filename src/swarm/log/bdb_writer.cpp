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
 *    \brief Defines an I/O interface for writing to files.
 *
 *
 */

#include "../common.hpp"

#include "../types/config.hpp"
#include "../plugin.hpp"

#include "log.hpp"
#include "writer.h"

#include <db_cxx.h>

namespace swarm {

using namespace gpulog;

/**
 * A writer that writes directly to Berkeley DB databases
 *
 *
 */
class bdb_writer : public swarm::writer
{
	// To do this, first we have to include the BDB in the CMake files
	// in process, we should make a gpulog::ilogstream out of the data
	// and read all the gpulog records and put them in the database
	// now the problem is what are we going to use for the key. Maybe we
	// can use a very simple key and later on we can do secondary databases
	// to make the indexes based on time and system id.
	Db db;
public:
	bdb_writer(const config& cfg):db(0,0){
		std::string fileName = cfg.require("log_output_db",std::string()) + ".db";
		db.open(NULL, fileName.c_str() , NULL, DB_RECNO, DB_CREATE | DB_TRUNCATE, 0);
	}
	virtual void process(const char *log_data, size_t length) {
		ilogstream stream(log_data,length);
		while(logrecord lr = stream.next()){
			Dbt key(NULL,0);
			Dbt data((void*)lr.ptr,lr.len());
			db.put(NULL,&key,&data,DB_APPEND);
		}
	}
	~bdb_writer(){
		db.close(0);
	}
};

writer_plugin_initializer< bdb_writer >
	bdb_writer_plugin("bdb", "This is the Berkeley DB writer");

}
