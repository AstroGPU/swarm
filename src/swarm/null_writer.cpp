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
#include "config.hpp"
#include "plugins.hpp"
#include "writer.h"

namespace swarm {

/**
 * null_writer plugin for use in
 * io.cpp
 *
 */
class null_writer : public swarm::writer
{
public:
	null_writer(const config& cfg){}
	virtual void process(const char *log_data, size_t length) {}
};

writer_plugin_initializer< null_writer >
	null_writer_plugin("null", "This is the dummy writer");

}
