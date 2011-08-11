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

#include "common.hpp"
#include <vector>
#include <algorithm> // for swap
#include <memory>
#include <iostream>
#include <fstream>
#include "writer.h"
#include "io.hpp"
#include "plugin_manager.hpp"

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

writer *writer::create(const config& cfg)
{
        std::auto_ptr<writer> w;

		std::string name = cfg.at("log writer");
		std::string plugin_name = "writer_" + name;

		try {
			w.reset( (writer*) (get_plugin(plugin_name)->create(cfg)) );
		}catch(plugin_not_found& e){
			ERROR("Log writer " + name + " not found.");
		}

        return w.release();
}
}
