/*************************************************************************
 * Copyright (C) 2010 by Young In Yeo and the Swarm-NG Development Team  *
 *                                                                       *
 * This program is free software; you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation; either version 3 of the License.        *
 *                                                                       *
 * This program is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 * GNU General Public License for more details.                          *
 *                                                                       *
 * You should have received a copy of the GNU General Public License     *
 * along with this program; if not, write to the                         *
 * Free Software Foundation, Inc.,                                       *
 * 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ************************************************************************/

/*! \file hermite_gpu_bpt.cpp
 * \brief Constructor and factory for gpu_hermite_bpt_integrator
 * 
*/

#include <cassert>
#include "swarm.h"
#include "hermite_gpu_bpt.h"
#include <cuda.h>
#include <cuda_runtime.h>

namespace swarm {

namespace hermite_gpu_bpt {

/*!
 * \brief Constructor for hermite gpu integrator
 *
 * Timestep, and blocksize are set from configuration class.
 * @param[in] cfg configuration class
 */
gpu_hermite_bpt_integrator::gpu_hermite_bpt_integrator(const config &cfg)
{
	if(!cfg.count("time step")) ERROR("Integrator gpu_hermite requires a timestep ('time step' keyword in the config file).");
	h = atof(cfg.at("time step").c_str());


	threadsPerBlock = cfg.count("threads per block") ? atoi(cfg.at("threads per block").c_str()) : 64;
}

/*!
 * \brief Factory to create double/single/mixed hermite gpu integrator based on precision
 *
 * @param[in] cfg configuration class
 *
 * @return        pointer to integrator cast to integrator*
*/
extern "C" integrator *create_gpu_hermite_bpt(const config &cfg)
{
	return new gpu_hermite_bpt_integrator(cfg);
}

} // end namespace hermite_gpu_bpt
} // end namespace swarm
