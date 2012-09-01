# Project Settings
PROJECT(Swarm_NG)
SET(Swarm_NG_VERSION_MAJOR 1)
SET(Swarm_NG_VERSION_MINOR 0)

# CMake Settings
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

# Directory Structure 
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET(swarm_INCLUDE_DIRS 
	${CMAKE_SOURCE_DIR}/src 
	${CMAKE_BINARY_DIR}/src
	)
SET(SAMPLES ${CMAKE_SOURCE_DIR}/samples) 
SET(TESTDIR ${CMAKE_SOURCE_DIR}/test) 



# All the cached variables and parameters that control
# the compilation go in this file
INCLUDE(variables.cmake)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)


# Find the libraries and packages
FIND_PACKAGE(CUDA REQUIRED)
FIND_PACKAGE(Boost REQUIRED COMPONENTS program_options regex)
FIND_PACKAGE(OpenMP)
FIND_PACKAGE(BDB) 

if(${CUDA_VERSION} VERSION_LESS ${REQUIRED_CUDA_VERSION})
	MESSAGE(SEND_ERROR "Your CUDA installation is outdated. 
		Swarm requires CUDA ${REQUIRED_CUDA_VERSION} or later.")
endif()

if(OPENMP_FOUND)
 set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
 set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
 set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


INCLUDE_DIRECTORIES(${swarm_INCLUDE_DIRS} ${Boost_INCLUDE_DIR})

 

# Macro to build other executables based on swarm library
MACRO(SWARM_ADD_EXECUTABLE target)
	ADD_EXECUTABLE(${target} ${ARGN})
	TARGET_LINK_LIBRARIES(${target} swarmng)
ENDMACRO(SWARM_ADD_EXECUTABLE target)

# Macro to build other libraries based on swarm library
MACRO(SWARM_ADD_LIBRARY target)
	ADD_LIBRARY(${target} ${ARGN})
	TARGET_LINK_LIBRARIES(${target} swarmng)
ENDMACRO(SWARM_ADD_LIBRARY target)

INCLUDE(doc.cmake)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/src/config.h)


# Include src
ADD_SUBDIRECTORY(src)


tutorial_python(TutorialPython py/tutorial.py)

INCLUDE(test.cmake)
INCLUDE(custom.cmake)

