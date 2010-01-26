#include "swarm.h"
#include "euler.h"
#include <cassert>

gpu_euler_integrator::gpu_euler_integrator(const config &cfg)
{
	if(!cfg.count("h")) ERROR("Integrator gpu_euler needs a timestep ('h' keyword in the config file).");
	h = atof(cfg.at("h").c_str());

	threadsPerBlock = cfg.count("threads per block") ? atoi(cfg.at("threads per block").c_str()) : 64;
}

// factory
extern "C" integrator *create_gpu_euler(const config &cfg)
{
	return new gpu_euler_integrator(cfg);
}
