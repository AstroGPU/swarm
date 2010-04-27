/*! \file hermite_gpu.h
 *  \brief declares gpu_hermite_integrator
 *
 *  Note that while this clas derivers from integrator, it does not use gpu_generic_integrator
*/

#ifndef integ_hermite_h__
#define integ_hermite_h__

#include "swarm.h"

namespace swarm {
namespace hermite_gpu {

/*!
 * \brief gpu_hermite_integrator class
 *
 * @tparam real_hi double
 * @tparam real_lo float for single and mixed, double for double
 */
template< typename real_hi, typename real_lo>
class gpu_hermite_integrator : public integrator
{
protected:
	//! time step 
	float h;
	//! precision 
	int prec;
	
	dim3 gridDim;
	//! blocksize
	int threadsPerBlock;

public:
	/*!
	 * \brief Constructor for hermite gpu integrator
	 *
	 * @param[in] cfg configuration class, read in from a file,  needs a timestep, precision, and block size.
	 */
	gpu_hermite_integrator(const config &cfg);

	/*!
	 * \brief host function to invoke a kernel (double precision) 
	 *
	 * @param[in,out] ens gpu_ensemble for data communication
	 * @param[in] dT destination time 
	 */
	void integrate(gpu_ensemble &ens, double dT);

};


} // end namespace hermite_gpu
} // end namespace swarm

#endif
