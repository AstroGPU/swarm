/*************************************************************************
 * Copyright (C) 2011 by Saleh Dindar and the Swarm-NG Development Team  *
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

/*! \file hermite_omp.hpp
 *   \brief Defines OpenMP implementation for PEC2 Hermite integrator.
 *
 */

#include "hermite_cpu.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif

namespace swarm { namespace cpu {

#ifdef _OPENMP
template< class Monitor >
class hermite_omp : public hermite_cpu<Monitor> {
	public:
	typedef hermite_cpu<Monitor> base;

	hermite_omp(const config& cfg): base(cfg){}
	virtual void launch_integrator() {
#pragma omp parallel for
		for(int i = 0; i < base::_ens.nsys(); i++){
			base::integrate_system(base::_ens[i]);
		}
	}


};
#endif

} } // Close namespaces
