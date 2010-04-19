#include "swarm.h"
#include "swarmio.h"
#include "swarmlog.h"
#include <memory>
#include <cmath>
#include <iostream>
#include <valarray>
#include <sstream>
#include <fstream>
#include <vector>

#define SWATCH_STOP(s)  { cudaThreadSynchronize(); (s).stop(); }
#define SWATCH_START(s) { (s).start(); }

int main(int argc, const char **argv)
{
	using namespace swarm;

	if(argc != 2)
	{
		std::cerr << "Usage: " << argv[0] << " <integrator.cfg>\n";
		return -1;
	}
	std::string icfgfn = argv[1];

	// load the ensemble
	cpu_ensemble ens;
	load_ensemble("data", ens);
	unsigned int nprint = std::min(2, ens.nsys());

//	gpulog::internal::dump_ttraits(body_set(ens, 1, 1, 3));
//	return -1;

	ens_writer out("output.bin");
	out << ens;

	// performance stopwatches
	stopwatch swatch_kernel, swatch_mem, swatch_temps, swatch_all;

	// set up the integrator
	SWATCH_START(swatch_all);
	config cfg;
	load_config(cfg, icfgfn);
	std::auto_ptr<integrator> integ(integrator::create(cfg));
	std::string runon = cfg.count("runon") ? cfg["runon"] : "gpu";
	bool ongpu;
	     if(runon == "gpu") { ongpu = true; }
	else if(runon == "cpu") { ongpu = false; }
	else { ERROR("The 'runon' configuration file parameter must be one of 'gpu' or 'cpu'"); }
	std::cerr << "Integrator: " << cfg["integrator"] << ", executing on the " << (ongpu ? "GPU" : "CPU") << "\n";

	// initialize output log
	std::string wcfg = cfg.count("output") ? cfg["output"] : "binary log.bin";
	std::cerr << "Output: " << wcfg << "\n";
	init_logs(wcfg);

	// set end times of integration, first output time, and snapshot interval
	double dT, Toutputstep;
	get_config(dT, cfg, "dT");
	get_config(Toutputstep, cfg, "output interval");
	for(int sys = 0; sys != ens.nsys(); sys++)
	{
		ens.time_end(sys) = ens.time(sys) + dT;
		ens.time_output(sys, 0) = ens.time(sys);	// output immediately on start
		ens.time_output(sys, 1) = Toutputstep;		// output interval
	}

	// log initialization
	hlog.alloc(1000000);
	gpulog::alloc_device_log("dlog", hlog.capacity());
//	hlog.attach_sink(w.get());

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
	flush_logs();

	// print out timings
	double us_per_sys_all = (swatch_all.getTime() / ens.nsys()) * 1000000;
	double us_per_sys_kernel = (swatch_kernel.getTime() / ens.nsys()) * 1000000;
	std::cerr << "Time per system (integration)   : " << us_per_sys_kernel << " us.\n";
	std::cerr << "Time per system (setup+integr.) : " << us_per_sys_all << " us.\n";
	std::cerr << "GPU/CPU memcpy time             : " << swatch_mem.getTime()*1000 << " ms.\n";
	std::cerr << "Internal state initialization   : " << swatch_temps.getTime()*1000 << " ms.\n";

	std::cerr << "Final time = (" << ens.time(0) << ", " << ens.time(1) << ",...)\n";
	// both the integrator & the ensembles are automatically deallocated on exit
	// so there's nothing special we have to do here.
	return 0;
}
