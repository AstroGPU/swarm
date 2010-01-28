#ifndef integ_mvs_h__
#define integ_mvs_h__

#include "swarm.h"

class gpu_mvs_integrator : public integrator
{
protected:
	float h;
	
	dim3 gridDim;
	int threadsPerBlock;

public:
	gpu_mvs_integrator(const config &cfg);

public:
	void integrate(gpu_ensemble &ens, double T, writer &w);
};

#endif
