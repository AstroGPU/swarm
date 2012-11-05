/*************************************************************************
 * Copyright (C) 2011 by Eric Ford and the Swarm-NG Development Team  *
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

/*! \file mvs_omp.hpp
 *   \brief Defines OpenMP implementation of mixed 
 *          variables symplectic propagator on CPU.
 *
 */

#include "mvs_cpu.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif

namespace swarm { namespace cpu {

#ifdef _OPENMP
/**
 *  *EXPERIMENTAL*: This class is not thoroughly tested.
 *  \ingroup experimental
 */
template< class Monitor >
class mvs_omp : public mvs_cpu<Monitor> {
	public:
	typedef mvs_cpu<Monitor> base;

	mvs_omp(const config& cfg): base(cfg){}
	virtual void launch_integrator() {
#pragma omp parallel for
		for(int i = 0; i < base::_ens.nsys(); i++){
			base::integrate_system(base::_ens[i]);
		}
	}


};
#endif

} } // Close namespaces
