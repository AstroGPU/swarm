INCLUDE_DIRECTORIES(.)

# Include all integrators
ADD_SUBDIRECTORY(integrators)
ADD_SUBDIRECTORY(propagators)


# Build swarm
CUDA_ADD_LIBRARY(swarm SHARED 
	swarm/plugins.cpp
	peyton/binarystream.cpp peyton/util.cpp 
	peyton/memorymap.cpp peyton/fakemmap.cpp
	swarm/snapshot.cpp swarm/integrator.cpp swarm/writer.cpp 
	swarm/utilities.cpp swarm/io.cpp   swarm/logmanager.cpp
	${SWARM_INTEGRATORS})

SET(COMPILE_FLAGS /MP)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
TARGET_LINK_LIBRARIES(swarm ${CUDA_CUDA_LIBRARY})


# Build the command line utilities
ADD_SUBDIRECTORY(utils)
ADD_SUBDIRECTORY(tutorials)
