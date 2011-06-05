/*************
 *  Author : Saleh Dindar
 *
 *
 */
#include <iostream>
#include "swarm.h"
#include <boost/program_options.hpp>
#include <boost/program_options/positional_options.hpp>
#include "utils.hpp"
#include "hp/hp.hpp"

using namespace swarm;
using namespace std;

int DEBUG_LEVEL  = 0;
const int LOGARITHMIC_BASE = 2;


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

	// Initialize ensemble on host to be used with GPU integration.
	DEBUG_OUTPUT(1,"Generate initial conditions and save it into ensemble");
	hp::defaultEnsemble reference_ensemble = generate_ensemble(cfg);

	DEBUG_OUTPUT(3, "Make a copy of ensemble" );
	hp::defaultEnsemble ens = reference_ensemble.clone() ; // Make a copy of the CPU ensemble for comparison


	// performance stopwatches
	stopwatch swatch_kernel_gpu, swatch_upload_gpu, swatch_download_gpu, swatch_temps_gpu, swatch_init_gpu;
	stopwatch swatch_all;

	swatch_all.start(); // Start timer for entire program
	srand(42u);    // Seed random number generator, so output is reproducible



	// Start GPU timers for initialization
	swatch_init_gpu.start();

	std::auto_ptr<hp::integrator> integ_gpu((hp::integrator*)integrator::create(cfg));
	cudaThreadSynchronize();  // Block until CUDA call completes
	swatch_init_gpu.stop();   // Stop timer for cpu initialization

	DEBUG_OUTPUT(1, "Upload ensemble to GPU" );
	cudaThreadSynchronize();   // Block until CUDA call completes
	swatch_upload_gpu.start(); // Start timer for copyg initial conditions to GPU
	integ_gpu->set_ensemble(ens);
	cudaThreadSynchronize();   // Block until CUDA call completes
	swatch_upload_gpu.stop();  // Stop timer for copyg initial conditions to GPU

	std::cout << "Time, Energy Conservation Error " << std::endl;

	for(double time = 0; time < duration ; ) {

		if((logarithmic > 1) && (time > 0)) interval = time * (logarithmic - 1);

		double step_size = min(interval, duration - time );
		integ_gpu->set_duration ( step_size  );

		DEBUG_OUTPUT(1, "Integrator ensemble on GPU" );
		swatch_kernel_gpu.start(); // Start timer for GPU integration kernel
		integ_gpu->integrate();  // Actually do the integration w/ GPU!			
		cudaThreadSynchronize();  // Block until CUDA call completes
		swatch_kernel_gpu.stop(); // Stop timer for GPU integration kernel  

		DEBUG_OUTPUT(2, "Check energy conservation" );
		double max_deltaE = find_max_energy_conservation_error(ens, reference_ensemble );

		time += step_size;

		std::cout << time << ", " << max_deltaE << std::endl;

	}

	/// CSV output for use in spreadsheet software 
	std::cout << "\n# Benchmarking times \n" 
		<< "Integration (ms), Integrator initialize (ms), "
		<< " Initialize (ms), GPU upload (ms), Download from GPU (ms) \n"
		<< swatch_kernel_gpu.getTime()*1000. << ",    "
		<< swatch_temps_gpu.getTime()*1000. << ", "
		<< swatch_init_gpu.getTime()*1000. << ", "
		<< swatch_upload_gpu.getTime()*1000. << ", "
		<< swatch_download_gpu.getTime()*1000. << ", "
		<< std::endl;

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

	config cfg;

	// Default configuration
	{
		cfg["integrator"] = "hp_hermite"; // Set to use a GPU integrator
		cfg["runon"]      = "gpu";         // Set to runon GPU
		cfg["time step"] = "0.0005";       // time step
		cfg["precision"] = "1";
		cfg["duration"] = "31.41592";
		cfg["nbod"] = "3";
		cfg["nsys"] = "960";
		cfg["blocksize"] = "128";
	}

	if(vm.count("cfg")){
		std::string icfgfn =  vm["cfg"].as<std::string>();
		load_config(cfg,icfgfn);
	}

	std::cout << "# Integrator:\t" << cfg["integrator"] << "\n"
		<< "# Time step\t" << cfg["time step"] << "\n"
		<< "# Min time step\t" << cfg["min time step"] << "\n"
		<< "# Max time step\t" << cfg["max time step"] << "\n"
		<< "# No. Systems\t" << cfg["nsys"] << "\n"
		<< "# No. Bodies\t" << cfg["nbod"] << "\n"
		<< "# Blocksize\t" << cfg["blocksize"] << "\n"
		<< std::endl;

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

	////////////// STABILITY TEST /////// 
	DEBUG_OUTPUT(1,"Initialize swarm library ");
	swarm::init(cfg);         // Initialize the library

	stability_test(cfg);


}
