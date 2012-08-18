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
#include "integrators/hermite_adap.hpp"
#include "monitors/stop_on_ejection.hpp"
#include "monitors/composites.hpp"
#include "swarm/gpu/gravitation_accjerk.hpp"
#include "monitors/log_transit.hpp"

typedef gpulog::device_log L;
using namespace swarm::monitors;
using namespace swarm::gpu::bppt;
using swarm::integrator_plugin_initializer;

integrator_plugin_initializer<
		hermite_adap< stop_on_ejection<L>  , GravitationAccJerk >
	> hermite_adap_plugin("hermite_adap");

integrator_plugin_initializer<
	        hermite_adap< stop_on_ejection_or_close_encounter<L>  , GravitationAccJerk > >
	hermite_adap_close_encounter_plugin("hermite_adap_close_encounter");


integrator_plugin_initializer<hermite_adap< log_transit<L>  , GravitationAccJerk > >
	hermite_adap_log_plugin("hermite_adap_transit");
