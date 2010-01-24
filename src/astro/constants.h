#ifndef _astro_constants_h
#define _astro_constants_h

namespace peyton {
/// Various physical and mathematical constants
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
// usability and backwards compatibility alias
namespace ctn = constants;
}

/// Inplace conversion of degrees to radians
#define RAD(x) { x = x*peyton::ctn::d2r; }
/// Inplace conversion of radians to degrees
#define DEG(x) { x = x/peyton::ctn::d2r; }

#define __peyton_constants peyton::constants

#endif
