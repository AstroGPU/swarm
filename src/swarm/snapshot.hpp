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

#include "ensemble.hpp"

#include <string>
#include <exception>
using std::exception;
using std::string;

namespace swarm {

class snapshot {
	private:
	snapshot(){}

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

	public:
	struct readfileexception : public exception {};
	struct writefileexception : public exception {};

	static defaultEnsemble load(const string& filename) throw (readfileexception);
	static void save(defaultEnsemble& ens, const string& filename)  throw (writefileexception);

};

}

