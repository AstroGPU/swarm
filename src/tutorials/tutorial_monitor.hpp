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

/** \file tutorial_monitor.hpp
 *  \brief A tutorial on how to make the monitor, for more info see @ref TutorialMonitor.
 *
 */

/*
 *  This is a simple tutorial used in doxygen pages
 *  should go through program2doxygen before it
 *  can be used by doxygen.
 *  
 *
 */
// @page TutorialMonitor Tutorial for making a monitor (stopper/logger)
// A “monitor” class is responsible for determining when the system state 
// should be logged and when the GPU should cease integrating a given system. 
// If the monitor determines that logging is needed, then the system’s current 
// state is written to a buffer in the GPU global device memory.
//
//
// Let's examine the tutorial included in this distribution @ref src/tutorials/tutorial_monitor.hpp.
// In this tutorial, we write a monitor that examines the distance of bodies to the origin.
// For planetary systems with barycenter or star at origin, this is generally means that
// the planet is is considered *ejected* and is no longer not considered part of the planetary system.
//
// If any planet is too far from origin, planet's information is written to a log and the system
// will be disabled. Disabling a
// system (by setting the system state to disabled) prevents the system from advancing
//  in time until the state is set back to active. Usually the user
// application should examine the system and decide if the system should be
// re-activated.
//
// Let's take a look at the code, first comes the preamble:

#pragma once
#include <limits>

const double POSITIVE_INFINITY = std::numeric_limits<double>::max();

namespace swarm { namespace monitors {

// Similar to a propagator, the parameters that need to be speficied
// in the configuration files should be stored in a separate struct because
// the main class is a template.
//
// For our simple monitor we only need to read in the maximum allowed distance 
// of planet to origin.
//
// Notice that the parameter is optional, if it is not provided, we use
// a default value of POSITIVE_INFINITY, meaning that by default there
// is no restriction.
struct TutorialMonitor_params 
{
	double max_distance_from_origin;
	
	TutorialMonitor_params(const config &cfg)
	{
		max_distance_from_origin = cfg.optional("max_distance_from_origin",POSITIVE_INFINITY);
	}
};

// Now the defining class of the monitor begins. The class should work with any implementation
// of logging system (There are two in Swarm-NG: host log and GPU log). So the log_t class is
// provided as a template parameter.
//
template<class log_t>
class TutorialMonitor 
{
	
// The parameter struct should be typedefed as params so the 
// generic integrator can find it.
public:
typedef TutorialMonitor_params params;

// These references are necessary for the monitor to function
// _params is our parameters, _sys is a reference to the system
// we are working on, and _log is a reference to the object
// the performs all the logging.
//
// These variable are provided in the constructor.
private:
params _params;
ensemble::SystemRef& _sys;
log_t& _log;

public:
GPUAPI TutorialMonitor(const params& p,ensemble::SystemRef& s,log_t& l)
	:_params(p),_sys(s),_log(l){}


// The monitor is implemented as a function object. It has a default
// function: to monitor the system. When it is called like a function
// with the system relative thread ID, it should perform the test
// and make changes as necessary. The parameter thread_in_system is
// provided for multi-threaded monitors. For single-threaded monitors
// We just use the thread 0 and don't do anything for other threads.
//
GPUAPI void operator () (const int thread_in_system) 
{ 
	if(thread_in_system == 0)
	{
		bool need_to_take_action = is_there_any_ejection();
		
		if(need_to_take_action)
		{
			log::system(_log, _sys);
			_sys.disable();
		}
	}
}
// The plan is simple, we test the system to see if there
// is any ejection (more details below). If there is any,
// we need to take action, otherwise our work is done
// We take two actions: 
//
// 1. Write current state of the system to the log.
// 2. Disable the system
// 
// The developer is free choose to do any of these two, or both.
// Here we show both actions.
//
//
// For the test, we write a loop over all the bodies and use
// the accessor method distance_to_origin to get the value.
// If the distance is greater than the provided parameter
// then we break out of the loop and report that there is
// an ejection.
// 
// Otherwise, the test returns false.
GPUAPI bool is_there_any_ejection()
{
	for(int b = 0; b < _sys.nbod(); b++){
		double distance_to_origin = _sys[b].distance_to_origin();
		if(distance_to_origin > _params.max_distance_from_origin)
			return true;
	}
	return false;
}

// This concludes the implementation of the monitor
// 
// The following three lines are just closing braces for class and namespace	
}; // end class TutorialMonitor
} } // end namespace monitor, end namespace swarm

