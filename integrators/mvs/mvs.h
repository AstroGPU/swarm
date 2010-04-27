/*! \file mvs.h
 * \brief declares gpu_mvs_integrator
 *
 * \todo implememnt gpu_mvs_integrator class
*/

#ifndef integ_mvs_h__
#define integ_mvs_h__

#include "swarm.h"

namespace swarm {

class gpu_mvs_integrator : public integrator
{
protected:
	float h;
	
	dim3 gridDim;
	int threadsPerBlock;

public:
	gpu_mvs_integrator(const config &cfg);

public:
	void integrate(gpu_ensemble &ens, double T);
};

} // end namespace swarm

#endif
