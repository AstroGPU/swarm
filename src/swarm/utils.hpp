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

/*! \file utils.hpp
 *   \brief Defines utility routines for common tasks that most users would need.
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
#include "ensemble_alloc.hpp"
#include "types/config.hpp"
#include <ostream>

/**
 * \file utils.hpp
 * \brief Utility routines that is used by applications
 *  
 * These are routines that are not part of the core swarm but
 * are used quite frequently in applications. For easier use
 * of swarm, these routines can be very useful.
 *
 * They range from input/output routines, generators and trivial statistical
 * data analysis.
 *
 */

#define $__(x,line) (std::cerr << __FUNCTION__ << ":" << line << ": " << #x <<  " = " << (x) << std::endl)
#define DEBUG_OUTPUT(level,message) ( (DEBUG_LEVEL >= level) ? (std::cerr << __FUNCTION__ << ":" << __LINE__ << ": " << message << std::endl) : std::cerr )

#define INFO_OUTPUT(level,message) ( (DEBUG_LEVEL >= level) ? (std::cerr <<  message ) : std::cerr )

/**
 * Pretty printing an ensemble::range_t, as a value with tolerance
 */
std::ostream& operator << (std::ostream& o, const swarm::ensemble::range_t& r);

/** 
 * Generate an ensemble of random generated star-centeric planetary systems.
 * The star is at origin and is stationary. The planets are on circular orbits.
 * The spacing factor is set to 1.4 by default but can be changed via 
 * cfg["spacing_factor"] config parameter.
 * Number of systems and number of bodies should be specified through
 * cfg["nsys"] and cfg["nbod"] parameters.
 *
 * \param cfg	Should containt "nsys" and "nbod", optionally "spacing_factor"
 */
swarm::hostEnsemble generate_ensemble(swarm::config& cfg);

/** 
 * Deprecated: Validate the configuration
 *
 * The criterias are very general and will cause problems.\todo This function
 * needs to be restructured or at least updated
 */
bool validate_configuration(swarm::config& cfg);

/**
 * Helper function to find the energy conservation error.
 * The total energy for each system in ens and reference_ens is
 * calculated individually. Then the difference in energy for each system 
 * is calculated and divided by energy in refrence_ens. The 
 * worst (highest) error rate is returned as the result.
 *
 * \param ens	 The ensemble to be examined for energy conservation
 * 	     (after integration)
 * \param reference_ens	  The ensemble that is used as a reference 
 * 		(initial conditions)
 */
double find_max_energy_conservation_error
	(swarm::ensemble& ens, swarm::ensemble& reference_ens ) ;


        swarm::ensemble::range_t energy_conservation_error_range(swarm::ensemble& ens, swarm::ensemble& reference_ensemble ) ;

/**
 * Pretty print selected values in a config data structure
 *
 */
void outputConfigSummary(std::ostream& o,swarm::config& cfg);

/**
 * Default configuration. Prepopulate a config data structure with
 * some pre-tested configuration values for non-optional parameters.
 */
swarm::config default_config() ;

/**
 * Compare ensembles and find their similarity.
 * This is used when two ensembles are supposed to be identical (e.g. computed
 * by the similar process). The difference in position and velocity and
 * time is calculated and returned in output variables. The function
 * returns false if the ensembles are not comparable (e.g. different number
 * of systems or bodies), otherwise the return value is true.
 *
 */
bool compare_ensembles( swarm::ensemble& e1, swarm::ensemble &e2 , double & pos_diff, double & vel_diff, double & time_diff ) ;

//! Basic helper function to count the number of systems with SYSTEM_DISABLED status 
int number_of_disabled_systems(swarm::defaultEnsemble ens) ;

