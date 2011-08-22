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

/*! \file swarm.h
 *   \brief Main header file for swarm library
 *
 *   Declares class swarm_error.  Methods declared here are defined in swarmlib.cpp
 *
 *   Declares swarm_error, ensemble, cpu_ensemble, gpu_ensemble, writer (abstract), 
 *   Also contains some memory management funcctions
*/
#pragma once

#include "common.hpp"
#include "types/ensemble.hpp"
#include "types/config.hpp"
#include "log/logmanager.hpp"
#include "integrator.hpp"
#include "plugin.hpp"
#include "utils.hpp"

/// The main namespace for the Swarm-NG library
namespace swarm {

/* Initialize the swarm library.  This function must be called before
 * any other */
inline void init(const config &cfg) { 
	swarm::log::manager::default_log()->init(cfg);
}


} // end namespace swarm
