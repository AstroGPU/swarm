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

/*! \file hermite_omp.cpp
 *   \brief Initializes the hermite OpenMP integrator plugin. 
 *
 */




#include "integrators/hermite_omp.hpp"
#include "monitors/log_time_interval.hpp"
#include "monitors/stop_on_ejection.hpp"
#include "monitors/composites.hpp"

//! Declare the host_log variable
typedef gpulog::host_log L;
using namespace swarm::monitors;
using swarm::integrator_plugin_initializer;
using namespace swarm::cpu;


#ifdef _OPENMP
integrator_plugin_initializer<
  hermite_omp< stop_on_ejection<L> >
	> hermite_omp_plugin("hermite_omp");
#endif

