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
 *  Text file format:
 *  \verbatim
 *   <nsys> <nbod>
 *   
 *   
 *   <time> <active>
 *        <mass> 
 *        <x> <y> <z>
 *        <vx> <vy> <vz>
 *
 *        <mass> 
 *        <x> <y> <z>
 *        <vx> <vy> <vz>
 *
 *        <mass> 
 *        <x> <y> <z>
 *        <vx> <vy> <vz>
 *
 *   <time> <active>
 *        <mass> 
 *        <x> <y> <z>
 *        <vx> <vy> <vz>
 *
 *        <mass> 
 *        <x> <y> <z>
 *        <vx> <vy> <vz>
 *
 *        <mass> 
 *        <x> <y> <z>
 *        <vx> <vy> <vz>
 *
 *     .
 *     .
 *     .
 *
 *  \endverbatim
 *
 *   The example above is for 3-body systems.
 *
 */
class snapshot {
	private:
	snapshot(){}

	public:
	struct header {
		int nsys, nbod;
	};
	struct sys {
		double time;
		double_int active;
	};
	struct body {
		double pos[3], vel[3], mass;
	};

	struct readfileexception : public std::exception {};
	struct writefileexception : public std::exception {};

	static defaultEnsemble load(const string& filename) throw (readfileexception);
	static defaultEnsemble load_text(const string& filename) throw (readfileexception);
	static void save(defaultEnsemble& ens, const string& filename)  throw (writefileexception);
	static void save_text(defaultEnsemble& ens, const string& filename)  throw (writefileexception);

};

}

