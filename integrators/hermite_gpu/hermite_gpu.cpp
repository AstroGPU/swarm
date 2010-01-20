#include "swarm.h"
#include "hermite_gpu.h"
#include <cassert>

gpu_hermite_integrator::gpu_hermite_integrator(const config &cfg)
{
	if(!cfg.count("h")) ERROR("Integrator gpu_hermite needs a timestep ('h' keyword in the config file).");
	h = atof(cfg.at("h").c_str());
}

// factory
extern "C" integrator *create_gpu_hermite(const config &cfg)
{
	return new gpu_hermite_integrator(cfg);
}
