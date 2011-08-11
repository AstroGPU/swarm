/*************
 *  Author : Saleh Dindar <saleh@cise.ufl.edu>, (c) 2011
 *
 * Copyright: See COPYING file that comes with this distribution
 *
 */
#include <iostream>
#include <boost/program_options.hpp>
#include <boost/program_options/positional_options.hpp>
#include "utils.hpp"
#include "swarm/integrator.hpp"
#include "swarm/stopwatch.h"

using namespace swarm;
using namespace std;

int DEBUG_LEVEL  = 0;
const int LOGARITHMIC_BASE = 2;

#define SYNC cudaThreadSynchronize()

void stability_test(config& cfg){
	if(!validate_configuration(cfg) ) {
		std::cerr << "Invalid configuration" << std::endl;
		return;
	}

	double duration = (cfg["duration"] != "") ? atof(cfg["duration"].c_str()) : 10 * M_PI ;        
	double interval = (cfg["interval"] != "") ? atof(cfg["interval"].c_str()) : (duration / 10 ) ; 
	double logarithmic = (cfg["logarithmic"] != "") ? atof(cfg["logarithmic"].c_str()) : 0 ; 

	if(interval < 1e-02 ) {
		std::cerr << "Interval is too small : " << interval << std::endl;
		return;
	}
	if(duration < interval ) {
		std::cerr << "Duration should be larger than interval : " << duration << std::endl;
		return;
	}

	DEBUG_OUTPUT(1,"Initialize swarm library ");
	swarm::init(cfg);

	DEBUG_OUTPUT(1,"Generate initial conditions and save it into ensemble");
	defaultEnsemble reference_ensemble = generate_ensemble(cfg);

	DEBUG_OUTPUT(2, "Make a copy of ensemble for energy conservation test" );
	defaultEnsemble ens = reference_ensemble.clone() ;

	// performance stopwatches
	stopwatch swatch_all;
	swatch_all.start(); // Start timer for entire program
	srand(42u);    // Seed random number generator, so output is reproducible

	std::auto_ptr<integrator> integ(integrator::create(cfg));
	SYNC;

	integ->set_ensemble(ens);
	SYNC;

	std::cout << "Time, Energy Conservation Error " << std::endl;

	for(double time = 0; time < duration ; ) {

		if((logarithmic > 1) && (time > 0)) interval = time * (logarithmic - 1);

		double step_size = min(interval, duration - time );
		integ->set_duration ( step_size  );

		DEBUG_OUTPUT(2, "Integrator ensemble on GPU" );
		integ->integrate();

		SYNC;
		DEBUG_OUTPUT(2, "Check energy conservation" );
		double max_deltaE = find_max_energy_conservation_error(ens, reference_ensemble );

		time += step_size;

		std::cout << time << ", " << max_deltaE << std::endl;

	}

	SYNC;
	swatch_all.stop();

	/// CSV output for use in spreadsheet software 
	std::cout << "\n# Integration CPU/GPU time: "
		<< swatch_all.getTime()*1000. << " ms " << std::endl;

}


int main(int argc,  char **argv)
{
	namespace po = boost::program_options;

	// Parse command line arguements (making it easy to compare)
	po::positional_options_description pos;
	po::options_description desc(std::string("Usage: ") + argv[0] + " \nOptions");

	desc.add_options()
		("interval,i", po::value<std::string>() , "Stability test intervals")
		("duration,d", po::value<std::string>() , "Duration of the integration")
		("logarithmic,l", po::value<std::string>() , "Produce times in logarithmic scale" )
		("help,h", "produce help message")
		("cfg,c", po::value<std::string>(), "Integrator configuration file")
		("verbose,v", po::value<int>(), "Verbosity level (debug output) ")
		;

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).
			options(desc).positional(pos).run(), vm);
	po::notify(vm);

	// Print help message if requested or invalid parameters
	if (vm.count("help")) { std::cout << desc << "\n"; return 1; }

	if (vm.count("verbose") ) DEBUG_LEVEL = vm["verbose"].as<int>();

	config cfg = default_config();

	if(vm.count("cfg")){
		std::string icfgfn =  vm["cfg"].as<std::string>();
		load_config(cfg,icfgfn);
	}

	outputConfigSummary(std::cout,cfg);

	if( ( cfg["blocksize"] != "" )  && (cfg["threads per block"] == "" )) 
		cfg["threads per block"] = cfg["blocksize"];

	if(vm.count("duration")) {
		cfg["duration"] = vm["duration"].as<std::string>();
	}
	if(vm.count("interval")) {
		cfg["interval"] = vm["interval"].as<std::string>();
	}

	if(vm.count("logarithmic")) {
		cfg["logarithmic"] = vm["logarithmic"].as<std::string>();
	}

	stability_test(cfg);
}
