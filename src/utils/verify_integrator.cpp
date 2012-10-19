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

/*! \file verify_integrator.cpp
 *   \brief Implements a utility to verify integrators against one another. 
 *
 */

#include <iostream>
#include "swarm/integrator.hpp"
#include <boost/program_options.hpp>
#include <boost/program_options/positional_options.hpp>
#include "utils.hpp"
#include "swarm/plugin_manager.hpp"
#include "swarm/snapshot.hpp"

const int DEBUG_LEVEL = 0;

#define SYNC cudaThreadSynchronize()

using namespace swarm;
using namespace std;


vector<config> parse_commandline_and_config(int argc, char* argv[]){
	namespace po = boost::program_options;
	po::positional_options_description pos;
	po::options_description desc(std::string("Tool to verify integrators against one another\nUsage: ") + argv[0] + " \nOptions");

	pos.add("config", -1);

	desc.add_options()
		("config,c", po::value<vector<string> >(), "Configuration for integrators to be compared")
		("basecfg,b", po::value<std::string>(), "Base config common to both integrators")
		("duration,d", po::value<std::string>() , "Duration of the integration")
		("input,i", po::value<std::string>(), "Input file")
		("output,o", po::value<std::string>(), "Output file")
		("help,h", "produce help message")
		("plugins,p", "list all of the plugins")
		;

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).
			options(desc).positional(pos).run(), vm);
	po::notify(vm);

	//// Respond to switches 
	//
	if (vm.count("help")) { std::cout << desc << "\n"; exit(1); }

	if (vm.count("plugins")) {
		cout << swarm::plugin_help_message;
		exit(0);
	}

	// Configurations loaded, first one is base
	vector<config> configs;

	if(vm.count("config") > 0){
		vector<string> config_fns = vm["config"].as< vector<string> >();
		configs.resize(config_fns.size()+1);
		for(int i = 0; i < config_fns.size() ; i++) 
			load_config(configs[i+1], config_fns[i]);
	}else{
		cerr << "No configuration specified " << endl;
		exit(-1);
	}

	if(vm.count("basecfg"))
		load_config(configs[0],vm["basecfg"].as<std::string>());
	else 
		configs[0] = default_config();

	const int cmd_to_config_len = 3;
	const char* cmd_to_config[cmd_to_config_len] = { "input", "output", "duration" };
	for(int i = 0; i < cmd_to_config_len; i++)
		if(vm.count(cmd_to_config[i]))
			configs[0][cmd_to_config[i]] = vm[cmd_to_config[i]].as<std::string>();
	
	
	return configs;
}

int main(int argc, char* argv[]){

	vector<config> cfgs = parse_commandline_and_config(argc,argv);
	config& basecfg = cfgs[0];
	cout << "Base config\n";
	outputConfigSummary(std::cout,basecfg);

	double duration = (basecfg["duration"] != "") ? atof(basecfg["duration"].c_str()) : 10 * M_PI ;        

	// Load/Generate the ensemble
	defaultEnsemble init_ens;
	if(basecfg["input"]!="") {
		cout << "Loading from " << basecfg["input"] << endl;
		init_ens = swarm::snapshot::load(basecfg["input"]);	
	}else{
		cout << "Generating new ensemble:  " << basecfg["nsys"] << ", " << basecfg["nbod"] << endl;
		init_ens = generate_ensemble(basecfg);
	}

	// Initialize Swarm
	swarm::init(basecfg);

	DEBUG_OUTPUT(0,"Reference configuration");
	// Integrate the base integrator
	defaultEnsemble ref = init_ens.clone();
	std::auto_ptr<integrator> integ(integrator::create(basecfg));
	integ->set_ensemble(ref);
	integ->set_duration ( duration  );
	integ->integrate();
	SYNC;
	double max_deltaE = find_max_energy_conservation_error(ref, init_ens );
	std::cout << "Max Energy Conservation Error =  " << max_deltaE << std::endl;
	// Save the ensemble
	if(basecfg["output"]!="") {
		cout << "Saving to " << basecfg["output"] << endl;
		swarm::snapshot::save(ref,basecfg["output"]);	
	}

	for(int i = 1; i < cfgs.size(); i++){
		DEBUG_OUTPUT(0,"Configuration case #" << i);
		outputConfigSummary(std::cout,cfgs[i]);
		defaultEnsemble ens = init_ens.clone();
		// Integrate the base integrator
		std::auto_ptr<integrator> integ(integrator::create(cfgs[i]));
		integ->set_ensemble(ens);
		integ->set_duration ( duration  );
		integ->integrate();
		SYNC;

		// TODO: what should be compared???
		/// Energy conservation error
		double max_deltaE = find_max_energy_conservation_error(ens, init_ens );
		std::cout << "Max Energy Conservation Error =  " << max_deltaE << std::endl;
		double max_diff = find_max_energy_conservation_error(ens, ref );
		std::cout << "Max Energy Difference with Reference =  " << max_diff << std::endl;
		if(cfgs[i]["output"]!="") {
			cout << "Saving to " << cfgs[i]["output"] << endl;
			swarm::snapshot::save(ens,cfgs[i]["output"]);	
		}
	}



	return 0;
}
