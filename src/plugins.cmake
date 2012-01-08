# CPU integrators
LIST(APPEND SWARM_PLUGINS integrators/hermite_cpu.cpp)

# GPU integrators
LIST(APPEND SWARM_PLUGINS integrators/hermite.cu)
LIST(APPEND SWARM_PLUGINS integrators/hermite_adap.cu)
LIST(APPEND SWARM_PLUGINS integrators/rkck.cu)

# Integrators based on generic integrator
LIST(APPEND SWARM_PLUGINS propagators/euler.cu)
#LIST(APPEND SWARM_PLUGINS propagators/mvs.cu)
#LIST(APPEND SWARM_PLUGINS propagators/verlet.cu)
#LIST(APPEND SWARM_PLUGINS propagators/midpoint.cu)

