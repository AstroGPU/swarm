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
	void integrate(gpu_ensemble &ens, double T);
};

#endif
