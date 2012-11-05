/*************************************************************************
 ************************************************************************
 * We would like to hear your feedback on using Swarm-NG libraries. 
 * Please send your comments on this tutorial and Swarm-NG to our mailing
 * list at: http://groups.google.com/group/swarm-ng.
 * 
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


// The generic is the generic integrator. It takes three parameters:
// an integrator, a monitor and an implementation of gravitational force
// calculation.
integrator_plugin_initializer< generic< TutorialPropagator, stop_on_ejection<L>, GravitationAccJerk > >
	tutorial_prop_plugin("tutorial"
			,"This is the integrator based on tutorial propagator");


// For complete listing of these two files look at \ref src/tutorials/tutorial_propagator.hpp and \ref src/tutorials/tutorial_propagator_plugin.cu
