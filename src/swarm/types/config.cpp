/*************************************************************************
 * Copyright (C) 2008-2010 by Mario Juric & Swarm-NG Development Team    *
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

/*! \file config.cpp
 *  \brief Implements several utility  functions for public interface for swarm libaray
 *
*/

#include "../common.hpp"
#include "../peyton/util.hpp"

#include "config.hpp"

using peyton::util::trim;


namespace swarm {

//! Initialize the config objects
void config::initialize(const char* const key_value_pairs[][2],const int n){
	$_(n);
	for(int i =0 ; i < n ; i++ )
		insert( std::pair<std::string,std::string>( key_value_pairs[i][0], key_value_pairs[i][1] ) ) ;
	
}

//! Load the configuration file and check for format errors
config config::load(const std::string &fn, config cfg)
{
	std::ifstream in(fn.c_str());
	if(!in) ERROR("Cannot open configuration file '" + fn + "'.");

	std::string line;
	int iline = 0;
	while(std::getline(in, line))
	{
		iline++;
		line = trim(line);
		if(line.empty()) { continue; }
		if(line[0] == '#') { continue; }

		size_t eqpos = line.find('=');
		if(eqpos == std::string::npos) ERROR("Error on line " + line + ": '=' sign expected.");

		std::string key = trim(line.substr(0, eqpos)), val = trim(line.substr(eqpos+1));

		cfg[key] = val;
	}
	return cfg;
}


}
