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

/*! \file hermite.cu
 *   \brief Initializes the hermite integrator plugins. 
 *
 */
#include "integrators/hermite_integrator.hpp"
#include "monitors/composites.hpp"
#include "monitors/stop_on_ejection.hpp"
#include "monitors/mercury_mce_stat.hpp"
#include "monitors/log_time_interval.hpp"
#include "swarm/gpu/gravitation_accjerk.hpp"
//#include "monitors/log_transit.hpp"
//#include "monitors/log_rvs.hpp"
#include "swarm/gpu/gravitation_acc.hpp"

//! Declare devide_log variable
typedef gpulog::device_log L;
using namespace swarm::monitors;
using namespace swarm::gpu::bppt;
using swarm::integrator_plugin_initializer;

//! Initialize the hermite integrator plugin for stop_on_ejection and gravitation acceleration
integrator_plugin_initializer<hermite< mce_stat<L> , GravitationAccJerk > > hermite_plugin("hermite");

//! Initialize the hermite integrator plugin for stop_on_ejection_or_close_encounter and gravitation acceleration
integrator_plugin_initializer<hermite< stop_on_ejection_or_close_encounter<L>  , GravitationAccJerk > >
	hermite_close_encounter_plugin("hermite_close_encounter");

//! Initialize the hermite integrator plugin for log_time_interval and gravitation acceleration
integrator_plugin_initializer<hermite< log_time_interval<L>  , GravitationAccJerk > >
	hermite_log_plugin("hermite_log");

//! Initialize the hermite integrator plugin for log_transition and gravitation acceleration
//integrator_plugin_initializer<hermite< log_transit<L>  , GravitationAccJerk > >
//	hermite_transit_plugin("hermite_transit");
// integrator_plugin_initializer<hermite< stop_on_ejection_or_close_encounter<L>  , GravitationAccJerk > >
// 	hermite_close_encounter_plugin("hermite_close_encounter");
// 
// //! Initialize the hermite integrator plugin for log_time_interval and gravitation acceleration
// integrator_plugin_initializer<hermite< log_time_interval<L>  , GravitationAccJerk > >
// 	hermite_log_plugin("hermite_log");
// 
// //! Initialize the hermite integrator plugin for log_transition and gravitation acceleration
// integrator_plugin_initializer<hermite< log_transit<L>  , GravitationAccJerk > >
// 	hermite_transit_plugin("hermite_transit");

#if __CUDA_ARCH__ >= 200
//integrator_plugin_initializer<hermite< log_rvs<L>  , GravitationAccJerk > >
//	hermite_rv_plugin("hermite_rv");
#endif

