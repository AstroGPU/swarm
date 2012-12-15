/*************************************************************************
 * Copyright (C) 2012 by Saleh Dindar and the Swarm-NG Development Team  *
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

/*! \file log2db.cpp
 *   \brief Implements a utility to that converts binary log files to
 *   databases and create indexes on top of them.
 *
 */


#include <boost/program_options.hpp>
#include <boost/program_options/positional_options.hpp>

#include <iostream>
#include <db_cxx.h>

#include "swarm/swarm.h"
#include "swarm/query.hpp"
#include "binary_reader.hpp"

using namespace swarm;
using namespace std;

namespace po = boost::program_options;
using boost::bind;
using boost::shared_ptr;
using gpulog::logrecord;

int recordsLimit = 1000;
int DEBUG_LEVEL = 0;

po::variables_map argvars_map;
string inputFileName, outputFileName;

void parse_commandline_and_config(int argc, char* argv[]){

	po::options_description desc("Usage:\n \tlog2db [options]\nOptions");
	desc.add_options()
			("input,i", po::value<std::string>(), "Binary log file from swarm binary_writer")
			("output,o", po::value<std::string>(), "Name of the database output file")
			("number,n", po::value<int>(), "Number of records to convert")
			("verbose,v", po::value<int>(), "Verbosity level")
			("dump,d", "Dump all the records up to the number")
			("quiet,q", "Suppress all messages")
			("help,h", "Help message")
			;



	po::variables_map &vm = argvars_map;
	po::store(po::command_line_parser(argc, argv).
			options(desc).run(), vm);
	po::notify(vm);

	//// Respond to switches
	//
	if (vm.count("help")) { std::cout << desc << "\n"; exit(1); }
	if (vm.count("verbose") ) DEBUG_LEVEL = vm["verbose"].as<int>();
	if (vm.count("quiet") ) DEBUG_LEVEL = -1;

	if(vm.count("input"))
		inputFileName = vm["input"].as<string>();
	else{
		cerr << "Name of input file is required" << endl;
		exit(2);
	}

	if(vm.count("output"))
		outputFileName = vm["output"].as<string>();
	else{
		cerr << "Name of output file is required" << endl;
		exit(2);
	}

	if(vm.count("number"))
		recordsLimit = vm["number"].as<int>();


}


template<typename T>
void put_in_dbt(const T& t, Dbt* data){
	data->set_flags(DB_DBT_APPMALLOC);
	data->set_size(sizeof(T));
	data->set_data(new T(t));
}


bool extract_from_ptr(void* ptr, size_t size, double& time, int& sys){
	logrecord l((char*)ptr);
	assert(l.len() == size);

	// Based on the implementation in query.cpp
	switch(l.msgid()){
	case 1: case 2: case 11: case 15: case 16:
		l >> time >> sys;
		return true;

	default:
		return false;
	}

}

/**
 * Extract the system ID from a log record
 *
 * @param secondary: The secondary parameter is the database handle for the secondary.
 * @param key      : The key parameter is a Dbt referencing the primary key.
 * @param data     : The data parameter is a Dbt referencing the primary data item.
 * @param result   : The result parameter is a zeroed Dbt in which the callback function
 *                   should fill in data and size fields that describe the secondary key or keys.
 */
int lr_extract_sysid(Db *secondary, const Dbt *key, const Dbt *data, Dbt *result) {
	double time = -1; int sys = -1;

	if(extract_from_ptr(data->get_data(), data->get_size(), time, sys)){
		put_in_dbt(sys, result);
		return 0;
	} else
		return DB_DONOTINDEX;
}

int lr_extract_evtid(Db *secondary, const Dbt *key, const Dbt *data, Dbt *result) {
	logrecord l((char*)data->get_data());
	assert(l.len() == data->get_size());

	put_in_dbt(l.msgid(), result);

	return 0;
}

int lr_extract_time(Db *secondary, const Dbt *key, const Dbt *data, Dbt *result) {
	double time = -1; int sys = -1;

	if(extract_from_ptr(data->get_data(), data->get_size(), time, sys)){
		put_in_dbt(time, result);
		return 0;
	} else
		return DB_DONOTINDEX;
}


/**
 * We need a strategy for reading for binary log files,
 * All the serialization functions that work on logrecords are optimized to work
 * with buffers. But we cannot allocate a buffer of 30GB, and we cannot partition it
 * either.
 *
 * The easiest way is to memory map a large portion of the file, pass it to the ilogstream and
 * let it run its course until there is not enough data. We catch the exception and memory map
 * the next section of the file and re-adjust the pointers and create ilogstream from scratch so
 * it points to the correct position to resume reading the logrecord. When we reach the end of the file
 * we have to make sure we set the pointers correctly for the ilogstream.
 *
 * But that is too complicated, we can make it even easier: instead of memory mapping, just read sections of file
 * into a fixed-size buffer and fix the pointers. Wrapping ilogstream seems to be too much trouble so we don't bother
 * with it and just write our on Java-style iterator: an object with hasNext and next methods.
 *
 * for testing it we can use a dump option that uses output_record function from the query.cpp and writes all the
 * records to the standard output.
 *
 * The easiest way to put the records into the file is to use time-system-evtid as the key and build some indices on
 * top of it.
 * However, the actual layout, indices and the primary keys depend on the data and our algorithms.
 */

using gpulog::logrecord;
using swarm::query::output_record;

/// main program
int main(int argc, char* argv[]){
	parse_commandline_and_config(argc,argv);

	// Open the binary file, we should be able to
	ifstream input(inputFileName.c_str(), ios::binary);
	binary_reader input_reader(input);
	if(!input_reader.validate())
		throw std::runtime_error("The input file is not valid");

	// Open the database files: primary, system_index, time_index, event_index
	// associate the indices with the primary

/*	DbEnv dbenv(0);
	dbenv.set_cachesize(0, 5 * 1024 * 1024 , 0);
	dbenv.set_data_dir(".");
	dbenv.open(".",  DB_CREATE | DB_INIT_LOG |
            DB_INIT_LOCK | DB_INIT_MPOOL |
            DB_INIT_TXN , 0);*/
	Db primary(NULL, 0), system_idx(NULL, 0), time_idx(NULL, 0), event_idx(NULL, 0);
	primary.open(NULL, (outputFileName+".p.db").c_str(), "primary", DB_BTREE, DB_CREATE, 0);
	system_idx.set_flags(DB_DUP | DB_DUPSORT);
	system_idx.open(NULL, (outputFileName+".sys.db").c_str(), "system_idx", DB_BTREE, DB_CREATE , 0);
	time_idx.set_flags(DB_DUP | DB_DUPSORT);
	time_idx.open(NULL, (outputFileName+".time.db").c_str(), "time_idx", DB_BTREE, DB_CREATE  , 0);
	event_idx.set_flags(DB_DUP | DB_DUPSORT);
	event_idx.open(NULL, (outputFileName+".evt.db").c_str(), "event_idx", DB_BTREE, DB_CREATE , 0);

	primary.associate(NULL, &system_idx,  &lr_extract_sysid, DB_IMMUTABLE_KEY);
	primary.associate(NULL, &time_idx  ,  &lr_extract_time , DB_IMMUTABLE_KEY);
	primary.associate(NULL, &event_idx ,  &lr_extract_evtid, DB_IMMUTABLE_KEY);



	for(int i=0; i < recordsLimit && input_reader.has_next(); i++) {
		// read one record from binary file
		logrecord l = input_reader.next();

		if(argvars_map.count("dump") > 0) {
			output_record(std::cout, l);
			std::cout << std::endl;
		}

		// insert it into the primary database, the secondary indices are automatically populated.
		Dbt key(&i,sizeof(int));
		Dbt data((void*)l.ptr,l.len());
		primary.put(NULL,&key,&data,0);

	}

	// Close all the databases
	primary.close(0);
	system_idx.close(0);
	time_idx.close(0);
	event_idx.close(0);
	//dbenv.close(0);

	// Close the binary file
}




