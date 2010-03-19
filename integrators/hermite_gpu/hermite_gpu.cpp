#include <cassert>
#include "swarm.h"
#include "hermite_gpu.h"
#include <cuda.h>
#include <cuda_runtime.h>

namespace swarm {

namespace hermite_gpu {


template<typename real_hi, typename real_lo>
gpu_hermite_integrator<real_hi,real_lo>::gpu_hermite_integrator(const config &cfg)
{
	if(!cfg.count("h")) ERROR("Integrator gpu_hermite needs a timestep ('h' keyword in the config file).");
	h = atof(cfg.at("h").c_str());

	if(!cfg.count("precision")) ERROR("Integrator gpu_hermite needs precision ('precision' keyword in the config file).");
	prec = atoi(cfg.at("precision").c_str());

	threadsPerBlock = cfg.count("threads per block") ? atoi(cfg.at("threads per block").c_str()) : 64;
}

// factory
extern "C" integrator *create_gpu_hermite(const config &cfg)
{
	if(!cfg.count("precision")) ERROR("Integrator gpu_hermite needs precision ('precision' keyword in the config file).");
	int prec = atoi(cfg.at("precision").c_str());
	assert((prec>=1)&&(prec<=3));
	if (prec==2)
	  return new gpu_hermite_integrator<float,float>(cfg);
	else if (prec==3)
	  return new gpu_hermite_integrator<double,float>(cfg);
	else // prec==1
	  return new gpu_hermite_integrator<double,double>(cfg);
}

} // end namespace hermite_gpu
} // end namespace swarm
