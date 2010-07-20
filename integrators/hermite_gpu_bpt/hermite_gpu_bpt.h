/*************************************************************************
 * Copyright (C) 2010 by Young In Yeo and the Swarm-NG Development Team  *
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

/*! \file hermite_gpu_bpt.h
 *  \brief declares gpu_hermite_bpt_integrator
 *
 *  Note that while this clas derivers from integrator, it does not use gpu_generic_integrator
*/

#ifndef integ_hermite_bpt_h__
#define integ_hermite_bpt_h__

#include "swarm.h"

namespace swarm {
namespace hermite_gpu_bpt {

/*!
 * \brief gpu_hermite_bpt_integrator class
 * computing only in double
 *
 */
class gpu_hermite_bpt_integrator : public integrator
{
protected:
	//! time step 
	float h;
	double dT;
	
	dim3 gridDim;
	//! blocksize
	int threadsPerBlock;
	//! size of required shared memory that should be dynamically allocated at launch time
	int shared_memory_size;
	dim3 threadDim;

public:
	/*!
	 * \brief Constructor for hermite gpu integrator
	 *
	 * @param[in] cfg configuration class, read in from a file,  needs a timestep, precision, and block size.
	 */
	gpu_hermite_bpt_integrator(const config &cfg);

	/*!
	 * \brief host function to invoke a kernel (double precision) 
	 *
	 * @param[in,out] ens gpu_ensemble for data communication
	 * @param[in] dT destination time 
	 */
	void integrate(gpu_ensemble &ens, double dT);

	template<int N>
	friend  struct kernel_launcher;

};


} // end namespace hermite_gpu_bpt
} // end namespace swarm

#endif
