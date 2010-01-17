#include "swarm.h"
#include "mvs.h"
#include <cassert>

gpu_mvs_integrator::gpu_mvs_integrator(const config &cfg)
{
	if(!cfg.count("h")) ERROR("Integrator gpu_mvs needs a timestep ('h' keyword in the config file).");
	h = atof(cfg.at("h").c_str());
}

// factory
extern "C" integrator *create_gpu_mvs(const config &cfg)
{
	return new gpu_mvs_integrator(cfg);
}
