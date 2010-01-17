#include "swarm.h"
#include "euler.h"

namespace gpu_euler_aux
{
	//
	// Wrap all aux. functions in a separate namespace, to avoid
	// collisions with equally named functions from other integrators.
	//

	__device__ float3 operator*(const float3 &a, const float &b)
	{
		return make_float3(a.x*b, a.y*b, a.z*b);
	}

	__device__ float3 operator+(const float3 &a, const float3 &b)
	{
		return make_float3(a.x+b.x, a.y+b.y, a.z+b.z);
	}

}

__constant__ ensemble gpu_euler_ens;
__global__ void gpu_euler_integrator_kernel(float dT, float h)
{
	using namespace gpu_euler_aux;

	ensemble &ens = gpu_euler_ens;
	int sys = threadId();
	if(sys >= ens.nsys()) { return; }

	float    T = ens.T(sys);
	float Tend = T + dT;

	// propagate the system until we match or exceed Tend
	while(T < Tend)
	{
		for(int i = 1; i != ens.nbod(); i++) // starting from 1, assuming 0 is the central body
		{
			// load
			double x = ens.x(sys, i);
			double y = ens.y(sys, i);
			double z = ens.z(sys, i);
			double vx = ens.vx(sys, i);
			double vy = ens.vy(sys, i);
			double vz = ens.vz(sys, i);
#if 0
			float E1 = (vx*vx + vy*vy + vz*vz) / 2. + 1. / sqrt(x*x + y*y + z*z);
#endif

			// compute acceleration
			float r2 = x*x + y*y + z*z;
			float aux = - 1.f / (r2 * sqrt(r2)) * h;
			float dvx = x * aux;
			float dvy = y * aux;
			float dvz = z * aux;

			// advance
			x += vx * h;
			y += vy * h;
			z += vz * h;
			vx += dvx;
			vy += dvy;
			vz += dvz;

			// store
			ens.x(sys, i) = x;
			ens.y(sys, i) = y;
			ens.z(sys, i) = z;
			ens.vx(sys, i) = vx;
			ens.vy(sys, i) = vy;
			ens.vz(sys, i) = vz;

#if 0
			float E2 = (vx*vx + vy*vy + vz*vz) / 2. + 1. / sqrt(x*x + y*y + z*z);
			float err = (E2-E1)/E1;
			printf("(%d, %d) dE/E = %f\n", sys, i, err);
#endif
		}

		T += h;
	}

	ens.T(sys) = T;
}

void gpu_euler_integrator::integrate(gpu_ensemble &ens, float dT)
{
	// Upload the kernel parameters
	if(ens.last_integrator() != this)
	{
		ens.set_last_integrator(this);
		configure_grid(gridDim, threadsPerBlock, ens.nsys());

		cudaMemcpyToSymbol(gpu_euler_ens, &ens, sizeof(gpu_euler_ens));
	}

	// execute the kernel
	gpu_euler_integrator_kernel<<<gridDim, threadsPerBlock>>>(dT, h);
}

