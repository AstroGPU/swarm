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

/*! \file snapshot.hpp
 *   \brief Implements load and save functions for ensemble files. 
 *
 */

#pragma once

#include "common.hpp"
#include "types/ensemble.hpp"

using std::string;

namespace swarm {

/**
 *  \page TextFormat Text Formats
 *
 *
 *  Text formats are input and output file formats used by load_text and save_text functions in \ref snapshot class. 
 *
 *  They are ASCII files containing textual representation of double precision floating point numbers.
 *  All quantities are double precision point numbers except for nsys, nbod, id and state, which are integers.
 *  State is defined in the ensemble:
 *    - 0 means active
 *    - 1 is disabled
 *    - other values are intermediate inactive states.
 *
 *  Text file format:
\verbatim
<number of systems> <number of bodies> <number of system attributes> <number of body attributes>


<system unique id> <system time> <system state>
<system attribute> <system attribute> ...
     <mass of the body> 
     <x> <y> <z>
     <velocity x> <velocity y> <velocity z>
	 <body attribute> <body attribute> ...

     <mass of the body> 
     <x> <y> <z>
     <velocity x> <velocity y> <velocity z>
	 <body attribute> <body attribute> ...

     <mass of the body> 
     <x> <y> <z>
     <velocity x> <velocity y> <velocity z>
	 <body attribute> <body attribute> ...


<system unique id> <system time> <system state>
<system attribute> <system attribute> ...
     <mass of the body> 
     <x> <y> <z>
     <velocity x> <velocity y> <velocity z>
	 <body attribute> <body attribute> ...

     <mass of the body> 
     <x> <y> <z>
     <velocity x> <velocity y> <velocity z>
	 <body attribute> <body attribute> ...

     <mass of the body> 
     <x> <y> <z>
     <velocity x> <velocity y> <velocity z>
	 <body attribute> <body attribute> ...


  .
  .
  .
 \endverbatim
 *
 *   The example above is for 3-body systems.
 *
 *   Explanation:
 *
 *   Header:
 *    - \<number of systems\>: Number of systems contained in this file. There should be exact number of blocks describing
 *    the properties of the system and bodies in the system after the header.
 *    - \<number of bodies\>: Number of bodies per system. All systems have equal number of bodies
 *    - \<number of system attributes\>: Number of general-purpose attributes per system. The purpose of the attributes is
 *    determined by the application
 *    - \<number of body attributes\> : Number of general-purpose attributes per body. The purpose of the attributes is
 *    determined by the application
 *
 *    System Block:
 *     - \<system unique id\>: Unique identifier for the system. The systems may be reduced or removed. The order of the systems
 *     may also change. This unique identifier is used to match systems between input/output files.
 *     - \<system time\> : Physical time of the system.
 *     - \<system state\> : Integration state of the system: 0 means active, 1 is disabled, other values are intermediate 
 *     inactive states that may be reactivated.
 *     - \<system attribute\>: A general-purpose attribute of the system. There are a specific number of these attributes and the order is important.
 *
 *    Body Block:
 *    - \<mass of the body\>: Mass of the body 
 *    - \<x\> \<y\> \<z\>: Three components of the position of the body
 *    - \<velocity x\> <velocity y\> \<velocity z\>: Three components of the velocity of the body
 *    - \<body attribute\>: a general-purpose attribute of the body, repeated on the same line and the order must be 
 *    preserved.
 *
 */

/*! Static class containing methods to load/save ensembles to file.
 *   
 *   This class supports load and save in text and binary formats.
 *   For binary format use:
 *     - ens = load(filename)
 *     - save(ens,filename)
 *
 *   For text format use:
 *     - ens = load_text(filename)
 *     - save_text(ens,filename)
 *
 *
 *  Binary file format: 
 *    HEADER | SYS | BODY | BODY | BODY | SYS | BODY | BODY | BODY | SYS | BODY | BODY | BODY | ...
 * 
 *   The example above is for an ensemble with 3 bodies per system. HEADER, SYS and BODY are structs
 *   defined inside the snapshot class.(refer to \ref header, \ref sys, \ref body).
 *   
 *   header containts number of systems and numebr of bodies. It also contains the number of attributes
 *   per system and number of attributes per body. Since the number of attributes are set at compile time
 *   the number of attributes in the file should be less than or equal to the number of the attributes that
 *   library supports.
 *
 *   For text format refere to \ref TextFormat
 *
 */
namespace snapshot {

	//! Data structure used in binary files. This is meant to be
	//! found at offset 0 of the file.
	struct header {
		int nsys, nbod, nsysattr, nbodattr;
	};
	//! Data structure used in binary files. parameters for a system.
	//! This comes right after the header and is followed by nbod number
	//! of body structs.
	struct sys {
		double time;
		int id;
		int state;
		double attribute[ensemble::NUM_SYS_ATTRIBUTES];
	};
	//! Data structure used in binary files. parameters for each body.
	//! nbod instances of this struct follows after 
	//! each sys data structure.
	struct body {
		double pos[3], vel[3], mass, attribute[ensemble::NUM_BODY_ATTRIBUTES];
	};

	//! Raised when an error encountered reading a text or binary file.
	//! Used in load and load_text
	struct readfileexception : public std::exception {
		string filename;
		int lineno;
		string message;
		readfileexception( const string& filename, const string& message, const int& lineno = 0): filename(filename),message(message),lineno(lineno) {}
		virtual ~readfileexception()throw(){}
		virtual const char* what() const throw() {
			return ("Error reading " + filename + " : " + message).c_str();
		}
	};

	//! Raised when an error encountered writing to a text or binary file.
	//! Used in save and save_text
	struct writefileexception : public std::exception {
		string filename;
		int lineno;
		string message;
		writefileexception( const string& filename, const string& message, const int& lineno = 0): filename(filename),message(message),lineno(lineno) {}
		virtual ~writefileexception()throw(){}
		virtual const char* what() const throw() {
			return ("Error writing " + filename + " : " + message).c_str();
		}
	};

	//! Load binary snapshot file
	defaultEnsemble load(const string& filename) 
		throw (readfileexception);

	//! Loads textual snapshot file
	defaultEnsemble load_text(const string& filename) 
		throw (readfileexception);

	//! Save the ensemble to a binary file
	void save(defaultEnsemble& ens, const string& filename) 
		throw (writefileexception);

	//! Save the ensemble as a text file
	void save_text(defaultEnsemble& ens, const string& filename) 
		throw (writefileexception);

}

}

