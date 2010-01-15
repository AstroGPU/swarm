#include "swarm.h"
#include "integrators.h"

////////// Utilities (should ultimately be moved elsewhere)

// Computes the global linear ID of the thread. Used from kernels.
// NOTE: Supports 3D grids with 1D blocks of threads
inline __device__ uint32_t threadId()
{
	const uint32_t id = ((blockIdx.z * gridDim.y + blockIdx.y) * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;
	return id;
}

////////// Euler integrator kernel

__constant__ ensemble gpu_euler_ens;
__global__ void gpu_euler_integrator_kernel(float dT, float h)
{
	ensemble &ens = gpu_euler_ens;
	int sys = threadId();
	if(sys >= ens.nsys()) { return; }
	
	float    T = ens.T(sys);
	float Tend = T + dT;

	// propagate the system until we match or exceed Tend
	while(T < Tend)
	{
		T += h;
	}

	ens.T(sys) = T;
}

void gpu_euler_integrator::integrate(gpu_ensemble &ens, float dT)
{
	// Upload the kernel parameters
	if(ens.last_integrator() != this)
	{
		configure_grid(gridDim, threadsPerBlock, ens.nsys());

		cudaMemcpyToSymbol(gpu_euler_ens, &ens, sizeof(gpu_euler_ens));
		ens.set_last_integrator(this);
	}

	// execute the kernel
	gpu_euler_integrator_kernel<<<gridDim, threadsPerBlock>>>(dT, h);
}
