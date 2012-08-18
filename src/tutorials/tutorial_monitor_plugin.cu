#include "tutorial_propagator.hpp"
#include "tutorial_monitor.hpp"
#include "swarm/gpu/gravitation_accjerk.hpp"

typedef gpulog::device_log L;
using namespace swarm::monitors;
using swarm::integrator_plugin_initializer;



integrator_plugin_initializer< generic< TutorialPropagator, TutorialMonitor<L>, GravitationAccJerk > >
	tutorial_prop_plugin("tutorial"
			,"This is the integrator based on tutorial propagator and tutorial monitor");


//
// For complete listing of these two files look at \ref src/tutorials/tutorial_monitor.hpp and \ref src/tutorials/tutorial_monitor_plugin.cu
