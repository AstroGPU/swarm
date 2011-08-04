LIST(APPEND SWARM_INTEGRATORS integrators/hermite_cpu.cpp)
LIST(APPEND SWARM_INTEGRATORS integrators/hermite.cu)
LIST(APPEND SWARM_INTEGRATORS integrators/rkck.cu)
LIST(APPEND SWARM_INTEGRATORS integrators/euler.cu)
LIST(APPEND SWARM_INTEGRATORS integrators/verlet.cu)
LIST(APPEND SWARM_INTEGRATORS integrators/midpoint.cu)

LIST(APPEND SWARM_INTEGRATORS propagators/hermite_propagator.cu)

