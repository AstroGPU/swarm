#include "swarm.h"

////////// Utilities

// Computes the global linear ID of the thread. Used from kernels.
// NOTE: Supports 3D grids with 1D blocks of threads
inline __device__ uint32_t threadId()
{
	const uint32_t id = ((blockIdx.z * gridDim.y + blockIdx.y) * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;
	return id;
}

