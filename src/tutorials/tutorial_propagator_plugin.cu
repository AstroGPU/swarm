/*************************************************************************
 * Copyright (C) 2009-2010 by Eric Ford & the Swarm-NG Development Team  *
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

/** \file tutorial_propagator_plugin.cu
 *  \brief A integrator plugin based on tutorial propagator. 
 *
 */


// This is the .cu file where you define the plugin
// 
//
#include "tutorial_propagator.hpp"
#include "monitors/stop_on_ejection.hpp"
#include "swarm/gpu/gravitation_accjerk.hpp"

typedef gpulog::device_log L;
using namespace swarm::monitors;
using swarm::integrator_plugin_initializer;



integrator_plugin_initializer< generic< TutorialPropagator, stop_on_ejection<L>, GravitationAccJerk > >
	tutorial_prop_plugin("tutorial"
			,"This is the integrator based on tutorial propagator");


// For complete listing of these two files look at \ref src/tutorials/tutorial_propagator.hpp and \ref src/tutorials/tutorial_propagator_plugin.cu
