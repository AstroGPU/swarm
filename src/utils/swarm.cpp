/*************************************************************************
 * Copyright (C) 2010 by Mario Juric      the Swarm-NG Development Team  *
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
/*! \file swarm.cpp
    Main program for integrating an ensemble of many N-body systems on a GPU

    Swarm assumes that N is small (3-10) and systems evolve Newtonian Gravity 
*/

#include "swarm/swarm.h" 
#include "swarm/log.hpp"
#include "utils.hpp"
#include <memory>
#include <iostream>

#define SWATCH_STOP(s)  { cudaThreadSynchronize(); (s).stop(); }
#define SWATCH_START(s) { (s).start(); }

int DEBUG_LEVEL = 0;

int main(int argc, const char **argv)
{
	if(argc != 2)
	{
		std::cerr << "Usage: " << argv[0] << " <integrator.cfg>\n";
		return -1;
	}
	std::string icfgfn = argv[1];

	// performance swatches
	stopwatch swatch_kernel, swatch_mem, swatch_temps, swatch_all;
	SWATCH_START(swatch_all);

	// load configuration
	swarm::config cfg;
	swarm::load_config(cfg, icfgfn);

	double duration = (cfg["duration"] != "") ? atof(cfg["duration"].c_str()) : 10 * M_PI ;        
	double interval = (cfg["interval"] != "") ? atof(cfg["interval"].c_str()) : (duration / 10 ) ; 
	double logarithmic = (cfg["logarithmic"] != "") ? atof(cfg["logarithmic"].c_str()) : 0 ; 

	if(interval < 1e-02 ) {
		std::cerr << "Interval is too small : " << interval << std::endl;
		return -1;
	}
	if(duration < interval ) {
		std::cerr << "Duration should be larger than interval : " << duration << std::endl;
		return -1;
	}


	std::string runon = cfg.count("runon") ? cfg["runon"] : "gpu";
	bool ongpu;
	if(runon == "gpu") { ongpu = true; }
	else if(runon == "cpu") { ongpu = false; }
	else { ERROR("The 'runon' configuration file parameter must be one of 'gpu' or 'cpu'");  return -1; }

//	outputConfigSummary(std::cout,cfg);

	// load the ensemble
	std::string ensprefix;
	swarm::get_config(ensprefix, cfg, "initial conditions");
	swarm::hostEnsemble ens;
	swarm::load_ensemble(ensprefix, ens);
	std::cout << "#Ensemble: \n#\tData Prefix: " << ensprefix 
		<< "\n#\tNo. Systems: " << ens.nsys() 
		<< "\n#\tNo. Bodies: " << ens.nbod() << std::endl;

	for(int i = 0; i < ens.nsys(); i++)
		ens.set_active(i);

	// initialize swarm -- this is required before calling any (non-utility) swarm library function
	swarm::init(cfg);

	// select the integrator to use
	std::auto_ptr<swarm::integrator> integ(swarm::integrator::create(cfg));

	DEBUG_OUTPUT(2,"Initializing integrator... ");
	SWATCH_START(swatch_temps);
	integ->set_default_log();
	integ->set_ensemble(ens);
	integ->set_duration(duration);

	// log initial conditions
	swarm::log::ensemble(hlog, ens);

	swarm::log::flush();

	SWATCH_STOP(swatch_temps);

	if(ongpu)
	{
		swarm::gpu::integrator* integ_gpu = (swarm::gpu::integrator*) integ.get();
		DEBUG_OUTPUT(2,"Uploading to GPU... ");
		SWATCH_START(swatch_mem);
		integ_gpu->upload_ensemble();
		SWATCH_STOP(swatch_mem);


		DEBUG_OUTPUT(2,"Integrating... ");
		SWATCH_START(swatch_kernel);
		integ_gpu->launch_integrator();
		SWATCH_STOP(swatch_kernel);

		DEBUG_OUTPUT(2,"Downloading data... ");
		SWATCH_START(swatch_mem);
		integ_gpu->download_ensemble();
		SWATCH_STOP(swatch_mem);
	}
	else
	{
		SWATCH_START(swatch_kernel);
		integ->integrate();				// integrate
		SWATCH_STOP(swatch_kernel);
	}
	for(int i =0; i < 20; i++){
		int sys = rand()%ens.nsys();
	}

	swarm::log::ensemble(hlog, ens);
	swarm::log::flush();
	SWATCH_STOP(swatch_all);
	DEBUG_OUTPUT(1,"Integration Complete");

	// print out timings
	double us_per_sys_all = (swatch_all.getTime() / ens.nsys()) * 1000000.0;
	double us_per_sys_kernel = (swatch_kernel.getTime() / ens.nsys()) * 1000000.0;
	std::cout << "# Time per system (integration)   : " << us_per_sys_kernel << " us.\n";
	std::cout << "# Time per system (setup+integr.) : " << us_per_sys_all << " us.\n";
	std::cout << "# GPU/CPU memcpy time             : " << swatch_mem.getTime()*1000.0 << " ms.\n";
	std::cout << "# Internal state initialization   : " << swatch_temps.getTime()*1000.0 << " ms.\n";

	return 0;
}
