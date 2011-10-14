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
#pragma once

#include "common.hpp"
#include "types/ensemble.hpp"

using std::string;

namespace swarm {

/**
 *  \page TextFormat
 *
 *
 *  Text formats are used by load_text and save_text functions in \ref snapshot class. 
 *
 *  The files are ASCII files containing textual representation of double precision floating point numbers.
 *  All quantities are double precision point numbers except for nsys, nbod, id and state, that are integers.
 *  State is defined in the ensemble. 0 means active, -1 is disabled, other values are intermediate 
 *  inactive states.
 *
 *  Text file format:
\verbatim
<nsys> <nbod>


<id> <time> <state>
     <mass> 
     <x> <y> <z>
     <vx> <vy> <vz>

     <mass> 
     <x> <y> <z>
     <vx> <vy> <vz>

     <mass> 
     <x> <y> <z>
     <vx> <vy> <vz>

<id> <time> <state>
     <mass> 
     <x> <y> <z>
     <vx> <vy> <vz>

     <mass> 
     <x> <y> <z>
     <vx> <vy> <vz>

     <mass> 
     <x> <y> <z>
     <vx> <vy> <vz>

  .
  .
  .
 \endverbatim
 *
 *   The example above is for 3-body systems.
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
 *   header containts number of systems and numebr of bodies.
 *
 *   For text format refere to \ref TextFormat
 *
 */
namespace snapshot {

	struct header {
		int nsys, nbod;
	};
	struct sys {
		double time;
		int id;
		int state;
	};
	struct body {
		double pos[3], vel[3], mass, attribute[NUM_PLANET_ATTRIBUTES];
	};

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

	/// Load binary snapshot file
	defaultEnsemble load(const string& filename) throw (readfileexception);
	/// Loads textual snapshot file
	defaultEnsemble load_text(const string& filename) throw (readfileexception);
	/// Save the ensemble to a binary file
	void save(defaultEnsemble& ens, const string& filename)  throw (writefileexception);
	/// Save the ensemble as a text file
	void save_text(defaultEnsemble& ens, const string& filename)  throw (writefileexception);

}

}

