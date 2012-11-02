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

/*! \file mvs.cu
 *   \brief Initializes the GPU version of the mixed variables symplectic propagator plugins.
 *
 */

#include "propagators/mvs.hpp"
#include "monitors/composites.hpp"
#include "monitors/stop_on_ejection.hpp"
#include "monitors/log_time_interval.hpp"
#include "swarm/gpu/gravitation_acc.hpp"

//! Declare device_log variable 
typedef gpulog::device_log L;
using namespace swarm::monitors;
using namespace swarm::gpu::bppt;
using swarm::integrator_plugin_initializer;

//! Initialize the integrator plugin for mvs propagator
integrator_plugin_initializer< generic< MVSPropagator, stop_on_ejection<L>, GravitationAcc > >
	mvs_prop_plugin("mvs"
			,"This is the integrator based on mvs propagator");

//! Initialize the integrator plugin for mvs propagator for close_encounter event
integrator_plugin_initializer< generic< MVSPropagator, stop_on_ejection_or_close_encounter<L>, GravitationAcc  > >
	mvs_prop_ce_plugin("mvs_close_encounter"
			,"This is the integrator based on mvs propagator, monitor stop_on_ejection_or_close_encounter");


