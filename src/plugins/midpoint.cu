/*************************************************************************
 * Copyright (C) 2011 by Saleh Dindar and the Swarm-NG Development Team  *
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

/*! \file midpoint.cu
 *   \brief Initializes the midpoint propagator integrator plugin. 
 *
 */

#include "propagators/midpoint.hpp"
#include "monitors/log_time_interval.hpp"
#include "swarm/gpu/gravitation_acc.hpp"

//! Declare device_log variable
typedef gpulog::device_log L;
using namespace swarm::monitors;
using namespace swarm::gpu::bppt;
using swarm::integrator_plugin_initializer;

//! Initialize the integrator plugin for midpoint propagator
integrator_plugin_initializer< generic< MidpointPropagator, stop_on_ejection<L>, GravitationAcc > >
	midpoint_prop_plugin("midpoint"
			,"This is the integrator based on midpoint propagator");

