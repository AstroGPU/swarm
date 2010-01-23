#include "swarm.h"
#include "swarmio.h"
#include <memory>
#include <cmath>
#include <iostream>
#include <valarray>
#include <sstream>
#include <fstream>
#include <vector>

#define SWATCH_STOP(s)  { cudaThreadSynchronize(); (s).stop(); }
#define SWATCH_START(s) { (s).start(); }

int main()
{
	// load the ensemble
	cpu_ensemble ens;
	load_ensemble("data", ens);
	unsigned int nprint = std::min(2, ens.nsys());

	ens_writer out("output.bin");
	out << ens;

	// performance stopwatches
	stopwatch swatch_kernel, swatch_mem, swatch_temps, swatch_all;

	// set up the integrator
	SWATCH_START(swatch_all);
	config cfg;
	load_config(cfg, "integrator.cfg");
	std::auto_ptr<integrator> integ(integrator::create(cfg));
	std::string runon = cfg.count("runon") ? cfg["runon"] : "gpu";
	bool ongpu;
	     if(runon == "gpu") { ongpu = true; }
	else if(runon == "cpu") { ongpu = false; }
	else { ERROR("The 'runon' configuration file parameter must be one of 'gpu' or 'cpu'"); }
	std::cerr << "Integrator: " << cfg["integrator"] << ", executing on the " << (ongpu ? "GPU" : "CPU") << "\n";
	// duration of integration
	double dT;
	get_config(dT, cfg, "dT");

	// perform the integration
	if(ongpu)
	{
		SWATCH_START(swatch_mem);
		gpu_ensemble gpu_ens(ens);				// upload to GPU
		SWATCH_STOP(swatch_mem);

		SWATCH_START(swatch_temps);
		integ->integrate(gpu_ens, 0.);				// initialize internal data structures
		SWATCH_STOP(swatch_temps);

		SWATCH_START(swatch_kernel);
		integ->integrate(gpu_ens, dT);				// integrate
		SWATCH_STOP(swatch_kernel);

		SWATCH_START(swatch_mem);
		ens.copy_from(gpu_ens);					// download to host
		SWATCH_STOP(swatch_mem);
	}
	else
	{
		SWATCH_START(swatch_temps);
		integ->integrate(ens, 0.);				// initialize internal data structures
		SWATCH_STOP(swatch_temps);

		SWATCH_START(swatch_kernel);
		integ->integrate(ens, dT);				// integrate
		SWATCH_STOP(swatch_kernel);
	}
	SWATCH_STOP(swatch_all);

	out << ens;

	// print out timings
	double us_per_sys_all = (swatch_all.getTime() / ens.nsys()) * 1000000;
	double us_per_sys_kernel = (swatch_kernel.getTime() / ens.nsys()) * 1000000;
	std::cerr << "Time per system (integration)   : " << us_per_sys_kernel << " us.\n";
	std::cerr << "Time per system (setup+integr.) : " << us_per_sys_all << " us.\n";
	std::cerr << "GPU/CPU memcpy time             : " << swatch_mem.getTime()*1000 << " ms.\n";
	std::cerr << "Internal state initialization   : " << swatch_temps.getTime()*1000 << " ms.\n";

	// both the integrator & the ensembles are automatically deallocated on exit
	// so there's nothing special we have to do here.
	return 0;
}
