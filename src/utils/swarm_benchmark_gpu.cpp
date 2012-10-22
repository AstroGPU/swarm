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

/*! \file swarm_benchmark_gpu.cpp
 *   \brief A utility for running performance benchmark on gpu
 *
 */

#include <iostream>
#include "swarm/integrator.hpp"
#include "swarm/stopwatch.h"
#include <boost/program_options.hpp>
#include <boost/program_options/positional_options.hpp>
#include "utils.hpp"
#include "swarm/logmanager.hpp"

using namespace swarm;
using namespace std;

int DEBUG_LEVEL  = 0;
#define SYNC cudaThreadSynchronize()

swarm::log::manager logman;

void run_integration(config& cfg, cpu_ensemble& reference_ensemble, const string& param, const string& value) {
	cfg[param] = value;
	if(!validate_configuration(cfg) ) {
		std::cerr << "Invalid configuration" << std::endl;
		return;
	}

  double duration = atof(cfg["duration"].c_str());

  // performance stopwatches
  stopwatch swatch_kernel_gpu, swatch_upload_gpu, swatch_download_gpu, swatch_temps_gpu, swatch_init_gpu;
  stopwatch swatch_all;

  swatch_all.start(); // Start timer for entire program
  srand(42u);    // Seed random number generator, so output is reproducible


  if( (param == "nsys") || (param == "nbod") ) {
	  DEBUG_OUTPUT(3, "Regenerating reference ensemble" );
	  reference_ensemble = generate_ensemble(cfg);
  }

  DEBUG_OUTPUT(3, "Make a copy of ensemble" );
  defaultEnsemble ens = reference_ensemble.clone() ;


  // Start GPU timers for initialization
  swatch_init_gpu.start();

  std::auto_ptr<gpu::integrator> integ_gpu((gpu::integrator*) integrator::create(cfg));
  integ_gpu->set_ensemble(ens);
  integ_gpu->set_duration ( duration );
  integ_gpu->set_log_manager( logman );
  SYNC;  // Block until CUDA call completes
  swatch_init_gpu.stop();   // Stop timer for cpu initialization

  DEBUG_OUTPUT(1, "Upload ensemble to GPU" );
  SYNC;   // Block until CUDA call completes
  swatch_upload_gpu.start(); // Start timer for copyg initial conditions to GPU
  integ_gpu->upload_ensemble();
  SYNC;   // Block until CUDA call completes
  swatch_upload_gpu.stop();  // Stop timer for copyg initial conditions to GPU

  DEBUG_OUTPUT(1, "Integrator ensemble on GPU" );

  swatch_kernel_gpu.start(); // Start timer for GPU integration kernel
  integ_gpu->launch_integrator();  // Actually do the integration w/ GPU!			
  SYNC;  // Block until CUDA call completes
  swatch_kernel_gpu.stop(); // Stop timer for GPU integration kernel  

  logman.flush();

  DEBUG_OUTPUT(1, "Download data to host" );
  swatch_download_gpu.start();  // Start timer for downloading data from GPU
  integ_gpu->download_ensemble();
  SYNC;      // Block until CUDA call completes
  swatch_download_gpu.stop();   // Stop timer for downloading data from GPU

  DEBUG_OUTPUT(2, "Check energy conservation" );
  double max_deltaE = find_max_energy_conservation_error(ens, reference_ensemble );

  /// CSV output for use in spreadsheet software 
  std::cout << param << ", ";
  std::cout << value << ",     ";
  std::cout << max_deltaE << ",    " ;
  std::cout << swatch_kernel_gpu.getTime()*1000. << ",    ";
  std::cout << swatch_temps_gpu.getTime()*1000. << ", ";
  std::cout << swatch_init_gpu.getTime()*1000. << ", ";
  std::cout << swatch_upload_gpu.getTime()*1000. << ", ";
  std::cout << swatch_download_gpu.getTime()*1000. << ", ";
  std::cout << std::endl;

}


void benchmark_loop(config& cfg, cpu_ensemble& ens, const string& param, const vector<string>& values) {

	for(vector<string>::const_iterator i = values.begin();  i != values.end(); i++){
		DEBUG_OUTPUT(1, "=========================================");
		run_integration(cfg, ens, param, *i);
	}

}

void benchmark_range(config& cfg, cpu_ensemble& ens, const string& param, int  from, int to , int increment ){

	for(int i = from; i <= to ; i+= increment ){
		std::ostringstream stream;
		stream <<  i;

		DEBUG_OUTPUT(1, "=========================================");
		run_integration(cfg, ens, param, stream.str() );
	}

}


int main(int argc,  char **argv)
{
	namespace po = boost::program_options;


	// Parse command line arguements (making it easy to compare)
	po::positional_options_description pos;
	po::options_description desc(std::string("Usage: ") + argv[0] + " \nOptions");

	pos.add("parameter", 1);
	pos.add("value", -1);


	desc.add_options()
		("value" , po::value<vector<string> >() , "Values to iterate over ")
		("parameter" , po::value<vector<string> >() , "Parameteres to benchmark ")
		("from", po::value<int>() , "from integer value")
		("to", po::value<int>() , "to integer value")
		("inc", po::value<int>() , "increment integer value")
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
		cfg["integrator"] = "hermite"; // Set to use a GPU integrator
		cfg["runon"]      = "gpu";         // Set to runon GPU
		cfg["time step"] = "0.0005";       // time step
		cfg["precision"] = "1";
		cfg["duration"] = "31.41592";
		cfg["nbod"] = "3";
		cfg["nsys"] = "960";
		cfg["blocksize"] = "16";
	}

	if(vm.count("cfg")){
		std::string icfgfn =  vm["cfg"].as<std::string>();
		load_config(cfg,icfgfn);
	}

	if( ( cfg["blocksize"] != "" )  && (cfg["threads per block"] == "" )) 
		cfg["threads per block"] = cfg["blocksize"];


	outputConfigSummary(std::cout,cfg);

	//////////////////////// BENCHMARKING /////////////////////// 

	// Initialize ensemble on host to be used with GPU integration.
	defaultEnsemble ens;

	DEBUG_OUTPUT(1,"Initialize swarm library ");
	swarm::init(cfg);         // Initialize the library

	DEBUG_OUTPUT(2,"Initialize logging system to null output ");
	{
		config logcfg; logcfg["log_writer"] = "null";
		logman.init(logcfg);
	}

	DEBUG_OUTPUT(1,"Generate initial conditions and save it into ensemble");
	ens = generate_ensemble(cfg);

	DEBUG_OUTPUT(1, "Column headers for CSV output ");
	std::cout << "Parameter, Value, Energy Conservation Error,  Integration (ms),";
	std::cout << " Integrator initialize (ms), Initialize (ms), GPU upload (ms), Download from GPU (ms)";
	std::cout << std::endl;

	if((vm.count("parameter") > 0) && (vm.count("value") > 0)) {

		string param = vm["parameter"].as< vector<string> >().front();
		vector<string> values = vm["value"].as< vector<string> >();

		benchmark_loop(cfg, ens, param,  values );

	} else if ((vm.count("parameter") > 0) &&(vm.count("from") > 0) && (vm.count("to") > 0)){

		int from = vm["from"].as<int>(), to = vm["to"].as<int>();
		int increment = vm.count("inc") > 0 ? vm["inc"].as<int>() : 1;
		string param = vm["parameter"].as< vector<string> >().front();

		benchmark_range(cfg, ens, param,  from, to , increment );

	} else {
		run_integration(cfg, ens, "" , "" );

	}
}

