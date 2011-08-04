INCLUDE_DIRECTORIES(.)


# Include all plugins
INCLUDE(plugins.cmake)

# Build swarm
CUDA_ADD_LIBRARY(swarm SHARED 
	swarm/plugin_manager.cpp
	peyton/binarystream.cpp peyton/util.cpp 
	peyton/memorymap.cpp peyton/fakemmap.cpp
	swarm/snapshot.cpp swarm/integrator.cpp 
	swarm/writer.cpp swarm/null_writer.cpp swarm/binary_writer.cpp 
	swarm/utilities.cpp swarm/io.cpp swarm/logmanager.cpp
	${SWARM_PLUGINS})
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})




# Build the command line utilities
ADD_SUBDIRECTORY(utils)
ADD_SUBDIRECTORY(tutorials)
