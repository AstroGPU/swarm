SET(GENERATE_FERMI TRUE CACHE BOOL "Whether to generate machine code for Fermi architecture (compute capability 2.0)")
SET(GENERATE_GT200 FALSE CACHE BOOL "Whether to generate machine code for GT200 architecture (compute capability 1.3)")
# Set CUDA Flags and options
SET(CUDA_NVCC_FLAGS 
	-Xcudafe --diag_suppress=subscript_out_of_range;
	-Xcudafe --diag_suppress=partial_override;
	-Xcudafe --diag_suppress=initialization_not_reachable;
	-w
	)
IF(${GENERATE_FERMI})
	SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} 
		-gencode arch=compute_20,code=sm_20;
		)
ENDIF()
IF(${GENERATE_GT200})
	SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} 
		-gencode arch=compute_13,code=sm_13;
		)
ENDIF()


SET(CUDA_TOOLKIT_ROOT_DIR /usr/local/cuda CACHE PATH "Path to CUDA SDK directory")

SET(NUM_PLANET_ATTRIBUTES 1 CACHE STRING "Number of attributes per planet [1..10]")
SET(NUM_SYSTEM_ATTRIBUTES 1 CACHE STRING "Number of attributes per system [1..10]")
SET(MIN_SHMEM_SIZE 17280 CACHE STRING "Minimum ammount of shared memory per block in bytes to be assumed at compile time")
SET(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/docs CACHE PATH "Where to put documentation output")




SET(MAX_NBODIES 6 CACHE STRING "Maximum number of bodies per system [3..10]")
# Automatically choosing the chunk size based on maximum number of bodies
# The constraint is that  max(nbod*3, nbod*(nbod-1)/2)*CHUNK_SIZE < 512
# the values used here are based on optimization on C2070

IF(${MAX_NBODIES} EQUAL 3)
	SET(OPTIMIZED_CHUNK_SIZE 16)
ELSEIF(${MAX_NBODIES} LESS 7)
	SET(OPTIMIZED_CHUNK_SIZE 8)
ELSEIF(${MAX_NBODIES} LESS 17)
	SET(OPTIMIZED_CHUNK_SIZE 4)
ELSEIF(${MAX_NBODIES} LESS 24)
	SET(OPTIMIZED_CHUNK_SIZE 2)
ELSEIF(${MAX_NBODIES} LESS 33)
	SET(OPTIMIZED_CHUNK_SIZE 1)
ELSE()
	MESSAGE(SEND_ERROR "Unsupported number of bodies: ${MAX_NBODIES}")	
ENDIF()
SET(ENSEMBLE_CHUNK_SIZE ${OPTIMIZED_CHUNK_SIZE} 
	CACHE STRING "Warpsize in ensemble for coalesced reads [1,4,8,16,32]")
SET(SHMEM_CHUNK_SIZE ${OPTIMIZED_CHUNK_SIZE} 
	CACHE STRING "Warpsize in shared memory for coalesced reads [1,4,8,16,32]")
MESSAGE("Optimal Chunk size:${OPTIMIZED_CHUNK_SIZE}")
MESSAGE("Current Shared memory Chunk size: ${SHMEM_CHUNK_SIZE}")
MESSAGE("Current Ensemble Chunk size: ${ENSEMBLE_CHUNK_SIZE}")

SET(CUDA_CUDA_LIBRARY /usr/lib/nvidia-current/libcuda.so CACHE PATH "Path to libcuda.so")
