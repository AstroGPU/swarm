#ifndef integ_hermite_adap_h__
#define integ_hermite_adap_h__

#include "swarm.h"

#ifndef LARGE_NUMBER
#define LARGE_NUMBER 1e32
#endif
#ifndef SMALL_NUMBER
#define SMALL_NUMBER 1e-32
#endif

namespace swarm {
namespace hermite_adap_gpu {

/**
 * gpu_hermite_adap_integrator class
 *
 * @tparam real_hi double
 * @tparam real_lo float for single and mixed, double for double
 */
template< typename real_hi, typename real_lo>
class gpu_hermite_adap_integrator : public integrator
{
protected:
	//! time step 
	real_hi h,stepfac;
	//! precision 
	int prec;
	
	dim3 gridDim;
	//! blocksize
	int threadsPerBlock;

public:
	/**
	 * Constructor for hermite gpu integrator
	 *
	 * @param[in] cfg configuration file needs a timestep, precision, stepfac, and block size.
	 */
	gpu_hermite_adap_integrator(const config &cfg);

	/**
	 * host function to invoke a kernel (double precision) 
	 *
	 * @param[in,out] ens gpu_ensemble for data communication
	 * @param[in] dT destination time 
	 */
	void integrate(gpu_ensemble &ens, double dT);

};


} // end namespace hermite_adap_gpu
} // end namespace swarm

#endif
