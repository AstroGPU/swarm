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

/*! \file text2snapshot.cpp
 *   \brief A utility for saving an ensemble. 
 *
 */

#include "swarm/snapshot.hpp"
#include <boost/program_options.hpp>
#include <boost/program_options/positional_options.hpp>
void parse_cmd(int argc, char* argv[], string& ifn, string& ofn){
	namespace po = boost::program_options;
	po::positional_options_description pos;
	po::options_description desc(std::string("Usage: ") + argv[0] + " \nOptions");

	desc.add_options()
		("input,i", po::value<std::string>(), "Input file")
		("output,o", po::value<std::string>(), "Output file")
		("help,h", "produce help message")
		;

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).
			options(desc).positional(pos).run(), vm);
	po::notify(vm);

	//// Respond to switches 
	//
	if ((vm.count("help")>0) || (vm.count("input") == 0) || (vm.count("output") == 0)) { std::cout << desc << "\n"; exit(1); }

	ifn = vm["input"].as<string>();
	ofn = vm["output"].as<string>();

}

int main(int argc,char * argv[]){
	string ifn, ofn;
	parse_cmd(argc,argv,ifn,ofn);
	swarm::defaultEnsemble ens = swarm::snapshot::load_text(ifn);
	swarm::snapshot::save(ens,ofn);
}

