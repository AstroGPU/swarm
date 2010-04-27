/*! \file mvs.cpp
 * \brief Constructor and factory for gpu_mvs_integrator
 * 
*/
#include "swarm.h"
#include "mvs.h"
#include <cassert>

namespace swarm {


/*!
 * \brief Constructor for mvs gpu integrator
 *
 * Timestep, precision, and blocksize are set from configuration class.
 * @param[in] cfg configuration class
 */
gpu_mvs_integrator::gpu_mvs_integrator(const config &cfg)
{
	if(!cfg.count("time step")) ERROR("Integrator gpu_mvs requires a timestep ('time step' keyword in the config file).");
	h = atof(cfg.at("time step").c_str());

	threadsPerBlock = cfg.count("threads per block") ? atoi(cfg.at("threads per block").c_str()) : 64;
}

/*!
 * \brief Factory to create mvs gpu integrator
 *
 * @param[in] cfg configuration class
 *
 * @return        pointer to integrator cast to integrator*
 */
extern "C" integrator *create_gpu_mvs(const config &cfg)
{
	return new gpu_mvs_integrator(cfg);
}
 
} // end namespace swarm
