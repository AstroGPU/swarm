SWARM_ADD_EXECUTABLE(tutorial_simple tutorial_simple.cpp)
SWARM_ADD_EXECUTABLE(tutorial_gpu tutorial_gpu.cpp)
SWARM_ADD_EXECUTABLE(montecarlo montecarlo.cpp kepler.cpp)
SWARM_ADD_EXECUTABLE(montecarlo_ecclimit montecarlo_ecclimit.cpp kepler.cpp)
SWARM_ADD_EXECUTABLE(montecarlo_mcmc_outputs montecarlo_mcmc_outputs.cpp kepler.cpp)



tutorial_file(TutorialSimple tutorial_simple.cpp )
tutorial_file(TutorialGPU tutorial_gpu.cpp    )
tutorial_file(TutorialPropagator tutorial_propagator.hpp tutorial_propagator_plugin.cu)
tutorial_file(TutorialMonitor tutorial_monitor.hpp tutorial_monitor_plugin.cu)
tutorial_file(TutorialIntegrator tutorial_integrator.hpp tutorial_integrator_plugin.cu)
