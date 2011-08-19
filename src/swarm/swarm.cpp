#include <iostream>
#include <boost/program_options.hpp>
#include <boost/program_options/positional_options.hpp>
#include "swarm/swarm.h"
#include "swarm/snapshot.hpp"
#include "stopwatch.h"

int DEBUG_LEVEL  = 0;
const int LOGARITHMIC_BASE = 2;
#define SYNC cudaThreadSynchronize()

using namespace swarm;
using namespace std;

void run_integration(config& cfg){
	if(!validate_configuration(cfg) ) {
		ERROR( "Invalid configuration" );
	}

	double duration = cfg.optional("duration", 10 * M_PI );
	double interval = cfg.optional("interval", duration / 10 ) ; 
	double logarithmic = cfg.optional("logarithmic", 0 ) ; 

	DEBUG_OUTPUT(1,"Initialize swarm library ");
	swarm::init(cfg);

	// Load/Generate the ensemble
	defaultEnsemble ref;
	if(cfg.valid("input") ) {
		DEBUG_OUTPUT(0,"Lading initial conditions from " << cfg["input"] );
		ref = swarm::snapshot::load(cfg["input"]);	
	}else{
		DEBUG_OUTPUT(0,"Generating new ensemble:  " << cfg["nsys"] << ", " << cfg["nbod"]);
		ref = generate_ensemble(cfg);
	}

	DEBUG_OUTPUT(2, "Make a copy of ensemble for energy conservation test" );
	defaultEnsemble ens = ref.clone();


	// Initialize Integrator
	std::auto_ptr<integrator> integ(integrator::create(cfg));
	integ->set_ensemble(ens);
	integ->set_duration ( duration  );
	SYNC;

	std::cout << "Time, Energy Conservation Error " << std::endl;

	// performance stopwatches
	stopwatch swatch_all;
	swatch_all.start(); // Start timer for entire program

	for(double time = 0; time < duration ; ) {

		if((logarithmic > 1) && (time > 0)) interval = time * (logarithmic - 1);

		double step_size = min(interval, duration - time );
		integ->set_duration ( step_size  );

		DEBUG_OUTPUT(2, "Integrator ensemble" );
		integ->integrate();

		SYNC;
		DEBUG_OUTPUT(2, "Check energy conservation" );
		double max_deltaE = find_max_energy_conservation_error(ens, ref );
		std::cout << time << ", " << max_deltaE << std::endl;

		time += step_size;

	}

	SYNC;
	swatch_all.stop();

	// Save the ensemble
	if(cfg.valid("output")) {
		cout << "Saving to " << cfg["output"] << endl;
		swarm::snapshot::save(ens,cfg["output"]);	
	}

	/// CSV output for use in spreadsheet software 
	std::cout << "\n# Integration time: "
		<< swatch_all.getTime()*1000. << " ms " << std::endl;

}

void parse_commandline_and_config(int argc, char* argv[], config& cfg){
	namespace po = boost::program_options;
	po::positional_options_description pos;
	po::options_description desc(std::string("Usage: ") + argv[0] + " \nOptions");

	desc.add_options()
		("duration,d", po::value<std::string>() , "Duration of the integration")
		("input,i", po::value<std::string>(), "Input file")
		("output,o", po::value<std::string>(), "Output file")
		("cfg,c", po::value<std::string>(), "Integrator configuration file")
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
	if (vm.count("verbose") ) DEBUG_LEVEL = vm["verbose"].as<int>();

	if (vm.count("plugins")) {
		cout << swarm::plugin::help;
		exit(0);
	}

	if(vm.count("cfg")){
		std::string icfgfn =  vm["cfg"].as<std::string>();
		cfg = config::load(icfgfn);
	}

	const int cmd_to_config_len = 3;
	const char* cmd_to_config[cmd_to_config_len] = { "input", "output", "duration" };
	for(int i = 0; i < cmd_to_config_len; i++)
		if(vm.count(cmd_to_config[i]))
			cfg[cmd_to_config[i]] = vm[cmd_to_config[i]].as<std::string>();

}

int main(int argc, char* argv[]){

	config cfg = default_config();
	parse_commandline_and_config(argc,argv,cfg);

	outputConfigSummary(std::cout,cfg);

	// Initialize Swarm
	swarm::init(cfg);

	run_integration(cfg);

	return 0;
}
