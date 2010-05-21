/*************************************************************************
 * Copyright (C) 2010 by Aaron Boley  and the Swarm-NG Development Team  *
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

/*! \file hermite_adap_gpu.h
 * \brief declares gpu_hermite_adap_integrator
 * 
 *  Note that while this class derivers from integrator, it does not use gpu_generic_integrator
*/

#ifndef integ_hermite_adap_h__
#define integ_hermite_adap_h__

#include "swarm.h"

namespace swarm {
namespace hermite_adap_gpu {

/**
 * gpu_hermite_adap_integrator class
 * \brief gpu_hermite_adap_integrator class
 *
 * @tparam real_hi double
 * @tparam real_lo float for single and mixed, double for double
 *
 * \todo Implement hermite_adap efficiently
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
