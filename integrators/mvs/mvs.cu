#include "swarm.h"
#include "mvs.h"

//
// FIXME: This is just a dummy placeholder (a copy of euler.cu).
//

namespace gpu_mvs_aux
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

__constant__ ensemble gpu_mvs_ens;
__global__ void gpu_mvs_integrator_kernel(double dT, float h)
{
	using namespace gpu_mvs_aux;

	ensemble &ens = gpu_mvs_ens;
	int sys = threadId();
	if(sys >= ens.nsys()) { return; }

	double    T = ens.time(sys);
	double Tend = T + dT;

	// propagate the system until we match or exceed Tend
	while(T < Tend)
	{
		for(int i = 0; i != ens.nbod(); i++)
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

	ens.time(sys) = T;
}

void gpu_mvs_integrator::integrate(gpu_ensemble &ens, double dT)
{
	ERROR("MVS integrator has not yet been implemented.");

	// Upload the kernel parameters
	if(ens.last_integrator() != this)
	{
		ens.set_last_integrator(this);
		configure_grid(gridDim, threadsPerBlock, ens.nsys());

		cudaMemcpyToSymbol(gpu_mvs_ens, &ens, sizeof(gpu_mvs_ens));
		if(dT == 0.) { return; }
	}

	// execute the kernel
	gpu_mvs_integrator_kernel<<<gridDim, threadsPerBlock>>>(dT, h);
}

