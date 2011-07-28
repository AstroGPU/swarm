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

/*! \file common.hpp
 *   \brief Common library headers between all files
 *
*/
#pragma once

#include <stdexcept>
#include <string>
#include <cstring>
#include <map>
#include <cassert>
#include <cmath>
#include <vector>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cux/cux.h>

#ifdef THROW_IS_ABORT
	#include <cassert>
	#include <cstring>
        #include <cstdio>
#endif

#ifndef __CUDACC__ // CUDA 2.2 C++ bug workaround
	#include <sstream>
        #include <valarray>
#endif

#define CUDADEVICETOUSE 1

#include <iostream>
#define $$(x) (std::cerr << __FILE__ << "(" << __FUNCTION__ << "):" << __LINE__ << " |> " << (x) << std::endl)
#define $_(x) (std::cerr << __FILE__ << "(" << __FUNCTION__ << "):" << __LINE__ <<  " " << (#x) << " = " << (x) << std::endl)
#define $$$ (std::cerr << __FILE__ << "(" << __FUNCTION__ << "):" << __LINE__ << " @@ " << std::endl)

#include "swarm/swarm_error.h"
