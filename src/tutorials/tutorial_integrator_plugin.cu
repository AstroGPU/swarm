// That was the definition of the plugin. Now we need to intantiate it as a 
// plugin so it is compiled and can be referenced in Swarm.
//
// For that we create a file with "cu" extension (for a GPU integrator). and
// include the definition of our integrator as well as plugin creating
// API.
//
// Here are the includes.
#include "tutorial_integrator.hpp"
#include "monitors/stop_on_ejection.hpp"
#include "swarm/gpu/gravitation_accjerk.hpp"

// Some aliasing and opening namespace can make our code more readable
typedef gpulog::device_log L;
using namespace swarm::monitors;
using swarm::integrator_plugin_initializer;

// First we have to define implementation for Monitor and Gravitational
// force calculation to our integrator to compile.
typedef TutorialIntegrator< stop_on_ejection<L>, GravitationAccJerk > TutorialIntegrator_Instantiated;

// Now that we have our integrator class, we use the Plug-in API to
// register our class with the name "tutorial_integrator" in the Swarm
// plugin system. You can use the name "tutorial_integrator" in the 
// config files to use the class e.g. integrator=tutorial_integrator
//
integrator_plugin_initializer< TutorialIntegrator_Instantiated >
	tutorial_integ_plugin("tutorial_integrator"
			,"This is the integrator based on tutorial integrator");

//
//
// For complete listing of these two files look at \ref src/tutorials/tutorial_integrator.hpp and \ref src/tutorials/tutorial_integrator_plugin.cu


/*************************************************************************
 ************************************************************************
 * We would like to hear your feedback on using Swarm-NG libraries. 
 * Please send your comments on this tutorial and Swarm-NG to our mailing
 * list at: http://groups.google.com/group/swarm-ng.
 * 
 ************************************************************************/
