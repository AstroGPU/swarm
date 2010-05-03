/*************************************************************************
 * Copyright (C) 2010 by Mario Juric  and the Swarm-NG Development Team  *
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
