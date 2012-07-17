/*************************************************************************
 * Copyright (C) 2010 by Saleh Dindar and the Swarm-NG Development Team  *
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
#pragma once

// TODO: Do we actually use this?
#include "../types/coalescedstructarray.hpp"
#include "bppt.hpp"

namespace swarm { namespace gpu { namespace bppt {

/**
 *  Unit type of the acceleration pairs shared array.
 *   
 *  This for each pair we keep an acc. The values are
 *  not final values, they are intermediate values calculated by
 *  calc_pair and should be accumulated using the correct algorithm
 *  to produce acceleration. CHUNK_SIZE can be 1. Usually it is
 *  set to 16 for optimizing coalesced reads from memory.
 *
 */
template<int W>
struct GravitationAccScalars {
	static const int CHUNK_SIZE = W;
	typedef double scalar_t; 

	double _acc[CHUNK_SIZE];

	// Accessors
	GENERIC double& acc() { return _acc[0];  }
};


/**
 *  Unit type of the acceleration and jerk pairs shared array.
 *   
 *  This for each pair we keep an acc and a jerk. The values are
 *  not final values, they are intermediate values calculated by
 *  calc_pair and should be accumulated using the correct algorithm
 *  to produce acceleration and jerk. CHUNK_SIZE can be 1. Usually it is
 *  set to 16 for optimizing coalesced reads from memory.
 *
 */
template<int W>
struct GravitationAccJerkScalars {
	static const int CHUNK_SIZE = W;
	typedef double scalar_t; 

	double _acc[CHUNK_SIZE];
	double _jerk[CHUNK_SIZE];

	// Accessors
	GENERIC double& acc() { return _acc[0];  }
	GENERIC double& jerk() { return _jerk[0];  }
};


/**
 * Helper function for calculating inner product
 */
GENERIC double inner_product(const double a[3],const double b[3]){
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

/// Helper function to convert an integer from 1..n*(n-1)/2 to a pair (first,second), this function returns the first element.
template<int nbod>
GENERIC int first ( int ij ){
	int i = nbod - 1 - ij / (nbod/2);
	int j = ij % (nbod/2);
	if (j < i) 
		return i;
	else 
		return nbod - 1 - i - nbod%2 + 1;
}

/// Helper function to convert an integer from 1..n*(n-1)/2 to a pair (first,second), this function returns the second element.
template<int nbod>
GENERIC int second ( int ij ){
	int i = nbod - 1 - ij / (nbod/2);
	int j = ij % (nbod/2);
	if (j < i) 
		return j;
	else 
		return nbod - 1 - j - nbod%2;
}


} } } // end of namespace bppt :: gpu :: swarm

