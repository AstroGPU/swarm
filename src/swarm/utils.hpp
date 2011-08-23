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

/*! \file utils.h
 *   \brief Utility routines for common tasks that most users would need.
 *
 *   Contains:
 *    - Generating ensembles
 *    - Comparing ensembles
 *    - Reporting functions
 *    - Some useful macros
 *
*/
#pragma once
#include "types/ensemble.hpp"
#include "types/config.hpp"
#include <ostream>

#define $__(x,line) (std::cerr << __FUNCTION__ << ":" << line << ": " << #x <<  " = " << (x) << std::endl)
#define DEBUG_OUTPUT(level,message) ( (DEBUG_LEVEL >= level) ? (std::cerr << __FUNCTION__ << ":" << __LINE__ << ": " << message << std::endl) : std::cerr )


std::ostream& operator << (std::ostream& o, const swarm::ensemble::range_t& r);

void generate_ensemble(swarm::config& cfg, swarm::cpu_ensemble& ens)  ;
bool validate_configuration(swarm::config& cfg);
double find_max_energy_conservation_error(swarm::ensemble& ens, swarm::ensemble& reference_ensemble ) ;
swarm::hostEnsemble generate_ensemble(swarm::config& cfg)  ;
void outputConfigSummary(std::ostream& o,swarm::config& cfg);
swarm::config default_config() ;
bool compare_ensembles( swarm::ensemble& e1, swarm::ensemble &e2 , double & pos_diff, double & vel_diff, double & time_diff ) ;
