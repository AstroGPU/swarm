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

/*! \file constants.hpp
 *   \brief Defines various constants for calculations. 
 *
 */

#ifndef _astro_constants_h
#define _astro_constants_h

namespace peyton {
//! Various physical and mathematical constants
namespace constants {

	const double pi = 3.14159265358979323846264338;
	const double pi2 = 2.*pi;
	const double twopi = 2.*pi;
	const double piby2 = pi/2.;
	const double halfpi = pi/2.;
	const double ln10 = 2.3025850929940456840;	// ln(10)

	const double d2r = pi/180.0;
	const double s2r = pi/(180.0*3600);

	// sunmass.h: Gaussian gravitational constant "k"
	const double gk = 0.01720209895; ///< k = sqrt(G) in units of Solar mass, days, AU
	const double gms = gk*gk;	///< k^2 = G in units of Solar mass, days, AU
	const double k2 = gms;		///< k^2 = G in units of Solar mass, days, AU

	// physical constants
	const double c = 299792458; ///< speed of light, [m/s]

	// unit conversion
	const double years = 3600*24*365; ///< number of seconds in a year
}
//! usability and backwards compatibility alias
namespace ctn = constants;
}

/// Inplace conversion of degrees to radians
#define RAD(x) { x = x*peyton::ctn::d2r; }
/// Inplace conversion of radians to degrees
#define DEG(x) { x = x/peyton::ctn::d2r; }

#define __peyton_constants peyton::constants

#endif
