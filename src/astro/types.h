/***************************************************************************
 *   Copyright (C) 2005 by Mario Juric   *
 *   mjuric@astro.Princeton.EDU   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef _astro_types_h
#define _astro_types_h

#include <cmath>

/**
	\file
	
	Declarations of libpeyton specific types.
	
	\todo Move sqr, sgn and sign templates to some other header in peyton::math namespace,
		that is, if they're needed at all
*/

namespace peyton {

typedef double MJD;
typedef double JD;
typedef double Degrees;
typedef double Radians;

template<typename T> inline T sqr(T x) { return x*x; }
template<typename T> inline T cube(const T &x) { return x*x*x; }
template<typename T> inline int sgn(T x) { return x > 0 ? 1 : (x < 0 ? -1 : 0); }
template<typename T> inline T sign(T a, T b) { return b >= 0 ? std::abs(a) : -std::abs(a); }

}

#define __peyton peyton

#endif
