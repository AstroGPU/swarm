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

/*! \file integrator.h
 *   \brief Declares integrator class 
 *
*/

#include "swarm_error.h"

/// The main namespace for the Swarm-NG library
namespace swarm {

typedef std::map<std::string, std::string> config;

/**
        \brief Abstract integrator interface

        Specific integrators must derive from this class (as well as define a
        create_XXXX(const config &cfg) factory function (where XXXX is the name
        of the integrator)). Instantiate them via a call to integrator::create(),
        where cfg["integrator"] is expected to contain the name of the integrator
        you wish to instantiate.
*/
class integrator
{
        public:
                virtual void integrate(gpu_ensemble &ens, double T)     //!< for GPU based integrators
                        { ERROR("Execution on GPU not supported by this implementation"); }
                virtual void integrate(cpu_ensemble &ens, double T)     //!< for CPU based integrators
                        { ERROR("Execution on CPU not supported by this implementation"); }

		// Destructor
                virtual ~integrator() {};       //!< has to be here to ensure the derived class' destructor is called (if it exists)

        public:
                /// Integrator factory functions (and supporting typedefs)
                static integrator *create(const config &cfg);

        protected:
                integrator() {};                /// hide the constructor.and force integrator instantiation with integrator::create
};

typedef integrator *(*integratorFactory_t)(const config &cfg);

}
