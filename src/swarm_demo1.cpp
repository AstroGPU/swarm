#include "swarm.h"
#include "swarmio.h"
#include <memory>
#include <cmath>
#include <iostream>
#include <valarray>
#include <sstream>
#include <fstream>
#include <vector>

int main(int argc, const char **argv)
{
        // Read configuration file
	if(argc != 2)
	{
		std::cerr << "Usage: " << argv[0] << " <integrator.cfg>\n";
		return -1;
	}
	config cfg;
	//	std::string icfgfn = argv[1];
	//	load_config(cfg, icfgfn);
	load_config(cfg, argv[1]);

	// load the initial conditions
	cpu_ensemble ens;
	load_ensemble("data", ens);

	// set up the integrator
	std::auto_ptr<integrator> integ(integrator::create(cfg));

	// Check that configuration file specifies which to use: CPU or GPU
	std::string runon = cfg.count("runon") ? cfg["runon"] : "gpu";
	bool ongpu;
	if(runon == "gpu") { ongpu = true; }
	else if(runon == "cpu") { ongpu = false; }
	else { ERROR("The 'runon' configuration file parameter must be one of 'gpu' or 'cpu'"); }
	std::cerr << "Integrator: " << cfg["integrator"] << ", executing on the " << (ongpu ? "GPU" : "CPU") << "\n";


	// set end times of integration, first output time, and snapshot interval
	double dT, Toutputstep;
	get_config(dT, cfg, "dT");
	//	get_config(Toutputstep, cfg, "output interval");
	for(int sys = 0; sys != ens.nsys(); sys++)
	{
	  ens.time_end(sys) = ens.time(sys) + dT;
	  //	  ens.time_output(sys, 0) = ens.time(sys);	// output immediately on start
	  //	  ens.time_output(sys, 1) = Toutputstep;		// output interval
	}

	// perform the integration
	if(ongpu)
	{
	  // upload to GPU
	  gpu_ensemble gpu_ens(ens);

	  // integrate
	  integ->integrate(gpu_ens, dT);				

	  // download to host
	  ens.copy_from(gpu_ens);					
	}
	else
	{
	  // integrate
	  integ->integrate(ens, dT);				
	}

	std::cerr << "Integration complete.\n";
	
	// both the integrator & the ensembles are automatically deallocated on exit
	// so there's nothing special we have to do here.
	return 0;
}



