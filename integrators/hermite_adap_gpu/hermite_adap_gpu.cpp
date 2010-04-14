#include "swarm.h"
#include "hermite_adap_gpu.h"
#include <cassert>

gpu_hermite_adap_integrator::gpu_hermite_adap_integrator(const config &cfg)
{
	if(!cfg.count("h")) ERROR("Integrator gpu_hermite needs a timestep ('h' keyword in the config file).");
	h = atof(cfg.at("h").c_str());

	if(!cfg.count("precision")) ERROR("Integrator gpu_hermite needs precision ('precision' keyword in the config file).");
	prec = atoi(cfg.at("precision").c_str());

	if(!cfg.count("stepfac")) ERROR("Integrator gpu_hermite needs a step factor ('stepfac' keyword in the config file).");
	stepfac = atof(cfg.at("stepfac").c_str());


	threadsPerBlock = cfg.count("threads per block") ? atoi(cfg.at("threads per block").c_str()) : 64;
}

// factory
extern "C" integrator *create_gpu_hermite_adap(const config &cfg)
{
	return new gpu_hermite_adap_integrator(cfg);
}
