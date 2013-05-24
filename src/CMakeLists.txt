INCLUDE_DIRECTORIES(.)

# Used for adding plugins to swarm
MACRO(ADD_PLUGIN mainfile id default_included text)
	SET(PLUGIN_${id} ${default_included} CACHE BOOL ${text})
	IF(${PLUGIN_${id}})
		LIST(APPEND SWARM_PLUGIN_FILES ${mainfile} ${ARGN})
		LIST(APPEND SWARM_PLUGINS ${id})

		SET(SWARM_PLUGIN_FILES ${SWARM_PLUGIN_FILES} PARENT_SCOPE)
		SET(SWARM_PLUGINS ${SWARM_PLUGINS} PARENT_SCOPE)
	ENDIF(${PLUGIN_${id}})
ENDMACRO()

# Include all plugins
INCLUDE(plugins.cmake)

MESSAGE("Enabled plugins are:")
FOREACH(p ${SWARM_PLUGINS})
	MESSAGE("\t${p}")
ENDFOREACH()
#MESSAGE("Compiled plugin files are: " ${SWARM_PLUGIN_FILES})


# Build swarm
INCLUDE_DIRECTORIES(${swarm_INCLUDE_DIRS} ${Boost_INCLUDE_DIR})
CUDA_ADD_LIBRARY(swarmng SHARED 
	swarm/plugin.cpp
	swarm/peyton/binarystream.cpp swarm/peyton/util.cpp 
	swarm/peyton/memorymap.cpp swarm/peyton/fakemmap.cpp
	swarm/snapshot.cpp swarm/integrator.cpp 
	swarm/log/writer.cpp swarm/log/null_writer.cpp 
    swarm/log/bdb_database.cpp
	swarm/log/io.cpp swarm/log/logmanager.cpp swarm/log/log.cpp
	swarm/gpu/device_settings.cpp
	swarm/types/config.cpp swarm/utils.cpp swarm/gpu/utilities.cu
	${SWARM_PLUGIN_FILES})
TARGET_LINK_LIBRARIES(swarmng ${Boost_LIBRARIES})
IF(BDB_FOUND)
	TARGET_LINK_LIBRARIES(swarmng ${BDB_LIBRARIES})
ENDIF(BDB_FOUND)

SWARM_ADD_EXECUTABLE(swarm swarm/swarm.cpp swarm/query.cpp swarm/bdb_query.cpp)




# Build the command line utilities
ADD_SUBDIRECTORY(utils)
ADD_SUBDIRECTORY(tutorials)
