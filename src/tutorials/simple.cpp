#include <iostream>
#include <boost/program_options.hpp>
#include <boost/program_options/positional_options.hpp>
#include "utils/utils.hpp"
#include "swarm/integrator.hpp"
#include "swarm/snapshot.hpp"
#include "swarm/logmanager.hpp"
#define SYNC cudaThreadSynchronize()

using namespace swarm;
using namespace std;

config default_config() {
	config cfg;
	cfg["integrator"] = "hermite"; // Set to use a GPU integrator
	cfg["time step"] = "0.0005";       // time step
	cfg["duration"] = "31.41592";
	cfg["nbod"] = "3";
	cfg["nsys"] = "960";
	cfg["blocksize"] = "16";
	return cfg;
}

void parse_commandline_and_config(int argc, char* argv[], config& cfg){
	namespace po = boost::program_options;
	po::positional_options_description pos;
	po::options_description desc(std::string("Usage: ") + argv[0] + " \nOptions");

	desc.add_options()
		("duration,d", po::value<std::string>() , "Duration of the integration")
		("input,i", po::value<std::string>(), "Input file")
		("output,o", po::value<std::string>(), "Output file")
		("help,h", "produce help message")
		("cfg,c", po::value<std::string>(), "Integrator configuration file")
		;

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).
			options(desc).positional(pos).run(), vm);
	po::notify(vm);

	//// Respond to switches 
	//
	if (vm.count("help")) { std::cout << desc << "\n"; exit(1); }

	if(vm.count("cfg")){
		std::string icfgfn =  vm["cfg"].as<std::string>();
		load_config(cfg,icfgfn);
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

	double duration = (cfg["duration"] != "") ? atof(cfg["duration"].c_str()) : 10 * M_PI ;        

	// Load/Generate the ensemble
	defaultEnsemble ref;
	if(cfg["input"]!="") {
		cout << "Loading from " << cfg["input"] << endl;
		ref = swarm::snapshot::load(cfg["input"]);	
	}else{
		cout << "Generating new ensemble:  " << cfg["nsys"] << ", " << cfg["nbod"] << endl;
		ref = generate_ensemble(cfg);
	}
	defaultEnsemble ens = ref.clone();

	// Initialize Swarm
	swarm::init(cfg);
	swarm::log::manager logman;
	logman.init(cfg["log output"]!="" ? cfg["log output"] : "null");

	// Initialize Integrator
	std::auto_ptr<integrator> integ(integrator::create(cfg));
	// TODO: detect gpu integrators
	{
		swarm::gpu::integrator* integ_gpu = (swarm::gpu::integrator*) integ.get();
		integ_gpu->set_log(logman.get_gpulog());
	}
	integ->set_ensemble(ens);
	integ->set_duration ( duration  );
	SYNC;

	// Integrate
	integ->integrate();
	logman.flush();
	SYNC;

	/// Energy conservation error
	double max_deltaE = find_max_energy_conservation_error(ens, ref );
	std::cout << "Max Energy Conservation Error =  " << max_deltaE << std::endl;

	// Save the ensemble
	if(cfg["output"]!="") {
		cout << "Saving to " << cfg["output"] << endl;
		swarm::snapshot::save(ens,cfg["output"]);	
	}

	return 0;
}
