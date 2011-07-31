INCLUDE_DIRECTORIES(.)

# Include all integrators
ADD_SUBDIRECTORY(integrators)
ADD_SUBDIRECTORY(propagators)


# Build swarm
CUDA_ADD_LIBRARY(swarm SHARED 
	${SWARM_INTEGRATORS}
	cux/cux.cpp 
	astro/BinaryStream.cpp astro/MemoryMap.cpp astro/Util.cpp 
	swarm/snapshot.cpp swarm/integrator.cpp swarm/writer.cpp 
	swarm/utilities.cpp swarm/io.cpp   swarm/logmanager.cpp)
TARGET_LINK_LIBRARIES(swarm ${CUDA_CUDA_LIBRARY})


# Build the command line utilities
ADD_SUBDIRECTORY(utils)
