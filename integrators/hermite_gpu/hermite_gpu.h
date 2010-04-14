#ifndef integ_hermite_h__
#define integ_hermite_h__

#include "swarm.h"

namespace swarm {
namespace hermite_gpu {


//class gpu_hermite_integrator : public integrator
template< typename real_hi, typename real_lo>
class gpu_hermite_integrator : public integrator
{
protected:
	float h;
	int prec;
	
	dim3 gridDim;
	int threadsPerBlock;

public:
	gpu_hermite_integrator(const config &cfg);

public:
	void integrate(gpu_ensemble &ens, double T);

};


} // end namespace hermite_gpu
} // end namespace swarm

#endif
