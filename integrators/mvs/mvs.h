/*************************************************************************
 * Copyright (C) 2010 by Mario Juric  and the Swarm-NG Development Team  *
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
