PROJECT(Swarm_NG)

# Basic Settings
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

SET(Swarm_NG_VERSION_MAJOR 0)
SET(Swarm_NG_VERSION_MINOR 2)


# Set CUDA Flags and options
SET(CUDA_NVCC_FLAGS 
	-gencode arch=compute_20,code=sm_20;
	-gencode arch=compute_13,code=sm_13;
	-Xcudafe --diag_suppress=subscript_out_of_range;
	-Xcudafe --diag_suppress=partial_override;
	-Xcudafe --diag_suppress=initialization_not_reachable;
	-w
	)

SET(CUDA_TOOLKIT_ROOT_DIR /usr/local/cuda CACHE PATH "Path to CUDA SDK directory")

SET(NUM_PLANET_ATTRIBUTES 1 CACHE STRING "Number of attributes per planet [1..10]")
SET(NUM_SYSTEM_ATTRIBUTES 1 CACHE STRING "Number of attributes per system [1..10]")
SET(MAX_NBODIES 6 CACHE STRING "Maximum number of bodies per system [3..10]")
SET(ENSEMBLE_CHUNK_SIZE 16 CACHE STRING "Warpsize in ensemble for coalesced reads [1,4,8,16,32]")
SET(SHMEM_CHUNK_SIZE 16 CACHE STRING "Warpsize in shared memory for coalesced reads [1,4,8,16,32]")
SET(MIN_SHMEM_SIZE 17280 CACHE STRING "Minimum ammount of shared memory per block in bytes to be assumed at compile time")
SET(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/docs CACHE PATH "Where to put documentation output")

# Find the libraries and packages
FIND_PACKAGE(CUDA REQUIRED)
FIND_PACKAGE(Boost REQUIRED COMPONENTS program_options regex)

SET(CUDA_CUDA_LIBRARY /usr/lib/nvidia-current/libcuda.so CACHE PATH "Path to libcuda.so")

# Macro to build other project based on swarm library
MACRO(SWARM_ADD_EXECUTABLE target)
	ADD_EXECUTABLE(${target} ${ARGN})
	INCLUDE_DIRECTORIES(${swarm_INCLUDE_DIRS} ${Boost_INCLUDE_DIR})
	TARGET_LINK_LIBRARIES(${target} swarmng ${Boost_LIBRARIES})
ENDMACRO(SWARM_ADD_EXECUTABLE target)

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
	add_custom_target(doc
		COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
		COMMENT "Generating API documentation with Doxygen" VERBATIM
		)
else(DOXYGEN_FOUND)
	add_custom_target(doc COMMENT "Doxygen not found, cannot make
	documentation")
endif(DOXYGEN_FOUND)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/src/config.h)

SET(swarm_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_BINARY_DIR}/src)

# Include src
ADD_SUBDIRECTORY(src)


SET(SAMPLES ${CMAKE_CURRENT_SOURCE_DIR}/samples) 
SET(TESTDIR ${CMAKE_CURRENT_SOURCE_DIR}/test) 

INCLUDE(test.cmake)
INCLUDE(custom.cmake)

