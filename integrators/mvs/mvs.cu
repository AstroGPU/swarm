/*************************************************************************
 * Copyright (C) 2010 by Mario Juric  and the Swarm-NG Development Team  *
 *                                                                       *
 * This program is free software; you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation; either version 3 of the License.        *
 *                                                                       *
 * This program is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 * GNU General Public License for more details.                          *
 *                                                                       *
 * You should have received a copy of the GNU General Public License     *
 * along with this program; if not, write to the                         *
 * Free Software Foundation, Inc.,                                       *
 * 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ************************************************************************/

/*! \file mvs.cu
 * \brief dummy placeholder until mvs is implemented
 *
 * \todo implement mvs GPU kernels (currently just placeholder)
 */

#include "swarm.h"
#include "mvs.h"

//
// FIXME: This is just a dummy placeholder (a copy of euler.cu).
//

namespace swarm {

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
#if 1
	ERROR("MVS integrator has not yet been implemented.");
#else
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
#endif
}

} // end namespace swarm
