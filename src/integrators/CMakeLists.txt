LIST(APPEND INTEGRATORS integrators/hermite_cpu.cpp)
LIST(APPEND INTEGRATORS integrators/hermite.cu)
LIST(APPEND INTEGRATORS integrators/rkck.cu)
LIST(APPEND INTEGRATORS integrators/euler.cu)
LIST(APPEND INTEGRATORS integrators/verlet.cu)
LIST(APPEND INTEGRATORS integrators/midpoint.cu)

# Give it to the parent
SET(SWARM_INTEGRATORS ${SWARM_INTEGRATORS} ${INTEGRATORS} PARENT_SCOPE)
