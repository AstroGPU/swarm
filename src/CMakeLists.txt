INCLUDE_DIRECTORIES(.)


# Include all plugins
INCLUDE(plugins.cmake)

# Build swarm
INCLUDE_DIRECTORIES(${swarm_INCLUDE_DIRS} ${Boost_INCLUDE_DIR})
CUDA_ADD_LIBRARY(swarmng SHARED 
	swarm/plugin.cpp
	swarm/peyton/binarystream.cpp swarm/peyton/util.cpp 
	swarm/peyton/memorymap.cpp swarm/peyton/fakemmap.cpp
	swarm/snapshot.cpp swarm/integrator.cpp 
	swarm/log/writer.cpp swarm/log/null_writer.cpp swarm/log/binary_writer.cpp 
	swarm/log/io.cpp swarm/log/logmanager.cpp
	swarm/gpu/device_settings.cpp
	swarm/types/config.cpp swarm/utils.cpp swarm/gpu/utilities.cu
	${SWARM_PLUGINS})

SWARM_ADD_EXECUTABLE(swarm swarm/swarm.cpp swarm/query.cpp)




# Build the command line utilities
#ADD_SUBDIRECTORY(utils)
ADD_SUBDIRECTORY(tutorials)
