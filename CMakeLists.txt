PROJECT(Swarm_NG)

# Basic Settings
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)


# Set CUDA Flags and options
SET(CUDA_VERBOSE_BUILD 0)
SET(CUDA_NVCC_FLAGS 
	-gencode arch=compute_20,code=sm_20;
	-gencode arch=compute_13,code=sm_13;
	-Xcudafe --diag_suppress=subscript_out_of_range;
	-Xcudafe --diag_suppress=partial_override;
	-Xcudafe --diag_suppress=initialization_not_reachable;
	-w
	)

include(CMakeUserConfig.txt)

SET(CUDA_TOOLKIT_ROOT_DIR /usr/local/cuda CACHE PATH "Path to CUDA SDK directory")

# Find the libraries and packages
FIND_PACKAGE(CUDA REQUIRED)
FIND_PACKAGE(Boost REQUIRED COMPONENTS program_options regex)

# Macro to build other project based on swarm library
MACRO(SWARM_ADD_EXECUTABLE target)
	ADD_EXECUTABLE(${target} ${ARGN})
	INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
	TARGET_LINK_LIBRARIES(${target} swarmng ${Boost_LIBRARIES})
ENDMACRO(SWARM_ADD_EXECUTABLE target)

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
	add_custom_target(doc
		${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
		COMMENT "Generating API documentation with Doxygen" VERBATIM
		)
endif(DOXYGEN_FOUND)

# Include src
ADD_SUBDIRECTORY(src)

