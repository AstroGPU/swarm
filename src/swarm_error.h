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

/*! \file swarm_error.h
 *   \brief Declares swarm error class
 *
*/
#include <stdexcept>
#include <string>

/// The main namespace for the Swarm-NG library
namespace swarm {

/**
        \brief Unrecoverable error exception.

        Throw an instance of this class to indicate an unrecoverable error
        was encountered. Do not throw it directly, but through the use of ERROR() macro.
*/
class swarm_error : public std::runtime_error
{
public:
        swarm_error(const std::string &msg) : std::runtime_error(msg) {}
        virtual ~swarm_error() throw() {};
};

#ifndef THROW_IS_ABORT
        #define ERROR(msg) throw swarm::swarm_error(msg);
#else
        #define ERROR(msg) { fprintf(stderr, "%s\n", std::string(msg).c_str()); abort(); }
#endif
}
