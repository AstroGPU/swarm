/*************************************************************************
 * Copyright (C) 2008-2010 by Mario Juric & Swarm-NG Development Team    *
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

/*! \file integrator.cpp
 *  \brief class integrator
 *
*/

#include <cuda_runtime_api.h>
#include <vector>
#include <algorithm> // for swap
#include <memory>
#include <iostream>
#include <dlfcn.h>
#include <fstream>
#include "swarm.h"
#include "swarmio.h"

namespace swarm {

/*!
   \brief Writer instantiation support

  @param[in] cfg configuration class
  @return writer
*/

/* Must use factory class to dynamically load integrator subclass
 * instead of using constructor. Done so that users can define their
 * own writers that the swarm code does not need to know about at
 * compile time
 */

writer *writer::create(const std::string &cfg)
{
        std::auto_ptr<writer> w;

        // try loading using a factory function
        void *me = dlopen(NULL, RTLD_LAZY);
        if(me == NULL)
        {
                ERROR(dlerror());
        }

        std::string name;
        std::istringstream ss(cfg);
        if(!(ss >> name))
                ERROR("Empty value for 'output' keyword in config file.");
        std::string factory_name = "create_writer_" + name;

        writerFactory_t factory = (writerFactory_t)dlsym(me, factory_name.c_str());
        if(factory)
        {
                std::string wcfg;
                getline(ss, wcfg);
                wcfg = trim(wcfg);
                w.reset(factory(wcfg));
        }
        else
        {
                ERROR("Writer " + name + " unknown (" + dlerror() + ").");
        }

        return w.release();
}
}
