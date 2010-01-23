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
