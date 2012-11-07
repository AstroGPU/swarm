// To use the monitor, we should write a plugin (or change any plugin to
// use our monitor.
//
// Here we examine @ref src/tutorials/tutorial_monitor_plugin.cu for an example:
//
// First we should include the source files that define implementation of
// 
// 1. A propagator that is used for the ODE integration
// 2. A monitor that logs and changes the state of the system as needed.
// 3. An algorithm for calculating gravitational forces.
//
#include "tutorial_propagator.hpp"
#include "tutorial_monitor.hpp"
#include "swarm/gpu/gravitation_accjerk.hpp"

// Some types need to be imported.
using gpulog::device_log;
using namespace swarm::monitors;
using swarm::integrator_plugin_initializer;


// We define the plugin in one line. The inner-working of the plugin system
// uses static initialization of object with constructors. We do not need
// to worry about that. We just use the integrator_plugin_initializer template
// as required. The template parameter to he integrator_plugin_initializer is
// Our integrator, composed of TutorialPropagator, TutorialMonitor<device_log> 
// (our monitor uses a GPU log object (device_log) since we are making a GPU integrator
// and GravitationAccJerk (TutorialPropagator uses jerk(third derivative of position
//  as well as acceleration).
//
// The two strings are just a name and a description. Choosing a name
// is important since the only way to access this plugin is through that
// name. 
integrator_plugin_initializer< generic< TutorialPropagator, TutorialMonitor<device_log>, GravitationAccJerk > >
	tutorial_prop_plugin("tutorial"
			,"This is the integrator based on tutorial propagator and tutorial monitor");
//
// To use this plugin, you should set the configuration parameter
// integrator=tutorial before creating the integrator. for example
// 
// @verbatim swarm --defaults integrate integrator=tutorial @endverbatim
//
//
//
// For complete listing of these two files look at \ref src/tutorials/tutorial_monitor.hpp and \ref src/tutorials/tutorial_monitor_plugin.cu
//
