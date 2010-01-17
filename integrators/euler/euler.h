#ifndef integ_euler_h__
#define integ_euler_h__

#include "swarm.h"

class gpu_euler_integrator : public integrator
{
protected:
	float h;
	
	dim3 gridDim;
	int threadsPerBlock;

public:
	gpu_euler_integrator(const config &cfg);

public:
	void integrate(gpu_ensemble &ens, float T);

	// No support for CPU execution. Note: we could make this function
	// transparently copy the CPU ensemble to the GPU, and back once
	// the integration is done
	void integrate(cpu_ensemble &ens, float T) { abort(); }
};

#endif
