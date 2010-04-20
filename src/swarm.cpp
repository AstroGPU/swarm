#include "swarm.h"
#include <memory>
#include <iostream>

#define SWATCH_STOP(s)  { cudaThreadSynchronize(); (s).stop(); }
#define SWATCH_START(s) { (s).start(); }

int main(int argc, const char **argv)
{
	if(argc != 2)
	{
		std::cerr << "Usage: " << argv[0] << " <integrator.cfg>\n";
		return -1;
	}
	std::string icfgfn = argv[1];
	std::cerr << "Configuration                   : " << icfgfn << "\n";

	// performance swatches
	stopwatch swatch_kernel, swatch_mem, swatch_temps, swatch_all;
	SWATCH_START(swatch_all);

	// load configuration
	swarm::config cfg;
	swarm::load_config(cfg, icfgfn);

	// load the ensemble
	std::string ensprefix;
	swarm::get_config(ensprefix, cfg, "ensemble");
	swarm::cpu_ensemble ens;
	swarm::load_ensemble(ensprefix, ens);
	std::cerr << "Ensemble                        : " << ensprefix << " (" << ens.nsys() << " systems, " << ens.nbod() << " bodies each).\n";

	// initialize swarm -- this is required before calling any (non-utility) swarm library function
	swarm::init(cfg);

	// select the integrator to use
	std::auto_ptr<swarm::integrator> integ(swarm::integrator::create(cfg));
	std::string runon = cfg.count("runon") ? cfg["runon"] : "gpu";
	bool ongpu;
	     if(runon == "gpu") { ongpu = true; }
	else if(runon == "cpu") { ongpu = false; }
	else { ERROR("The 'runon' configuration file parameter must be one of 'gpu' or 'cpu'"); }
	std::cerr << "Integrator                      : " << cfg["integrator"] << ", executing on the " << (ongpu ? "GPU" : "CPU") << "\n";

	// set end times of integration, first output time, and snapshot interval
	double dT, Toutputstep;
	swarm::get_config(dT, cfg, "dT");
	swarm::get_config(Toutputstep, cfg, "output interval");
	for(int sys = 0; sys != ens.nsys(); sys++)
	{
		ens.time_end(sys) = ens.time(sys) + dT;
		ens.time_output(sys, 0) = ens.time(sys);	// output immediately on start
		ens.time_output(sys, 1) = Toutputstep;		// output interval
	}

	// perform the integration
	std::cerr << "\nIntegrating. ";
	if(ongpu)
	{
		SWATCH_START(swatch_mem);
		swarm::gpu_ensemble gpu_ens(ens);			// upload to GPU
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
	std::cerr << "Done.\n\n";

	// print out timings
	double us_per_sys_all = (swatch_all.getTime() / ens.nsys()) * 1000000;
	double us_per_sys_kernel = (swatch_kernel.getTime() / ens.nsys()) * 1000000;
	std::cerr << "Time per system (integration)   : " << us_per_sys_kernel << " us.\n";
	std::cerr << "Time per system (setup+integr.) : " << us_per_sys_all << " us.\n";
	std::cerr << "GPU/CPU memcpy time             : " << swatch_mem.getTime()*1000 << " ms.\n";
	std::cerr << "Internal state initialization   : " << swatch_temps.getTime()*1000 << " ms.\n";

	return 0;
}
