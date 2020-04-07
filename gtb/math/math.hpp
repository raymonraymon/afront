
/*

Copyright 2007 University of Utah


This file is part of Afront.

Afront is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published
by the Free Software Foundation; either version 2 of the License,
or (at your option) any later version.

Afront is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA

*/


#ifndef GTB_MATH_INCLUDED
#define GTB_MATH_INCLUDED

#include <gtb/common.hpp>
#include <gtb/math/matrix4.hpp>
#include <gtb/math/mathtools.h>
#include <gtb/math/mynr.h>
#include <gtb/math/arraymat.h>
#include <gtb/math/functions.h>
#include <gtb/math/epsilon.hpp>
#include <gtb/graphics/vector3.hpp>
#include <gtb/math/poly.h>

#ifndef WIN32
//#	include <stdlib.h>
#	ifdef __GNUC__
#		define __USE_ISOC99 1
#	endif // __GNUC__
#	include <math.h>
#	include <algorithm>
#endif // WIN32

#ifndef isfinite
#	if defined(WIN32)
#		define isfinite _finite
#	elif defined(__sgi)
#		define isfinite finite
#       elif defined(__APPLE__)
#               define isfinite finite
#	else
#		define isfinite __finite
#	endif // WIN32
#endif	// isfinite


GTB_BEGIN_NAMESPACE


#ifndef M_PI
#define M_E             2.7182818284590452354
#define M_LOG2E         1.4426950408889634074
#define M_LOG10E        0.43429448190325182765
#define M_LN2           0.69314718055994530942
#define M_LN10          2.30258509299404568402
#define M_PI            3.14159265358979323846
#define M_PI_2          1.57079632679489661923
#define M_PI_4          0.78539816339744830962
#define M_1_PI          0.31830988618379067154
#define M_2_PI          0.63661977236758134308
#define M_2_SQRTPI      1.12837916709551257390
#define M_SQRT2         1.41421356237309504880
#define M_SQRT1_2       0.70710678118654752440
#define M_SQRT3         1.73205080756887729352745
#endif // M_PI


#define M_PI_OVER_180  0.01745329251994329576
#define M_180_OVER_PI 57.29577951308232087684


#define POINTS_PER_INCH 72
#define CM_PER_INCH 2.54
#define MM_PER_INCH 25.4


double deg_to_rad(double x);
double rad_to_deg(double x);

/*
 * convert spherical to euclidian coordinates and back
 */
/*typedef tVector3<float> Vector3f;
typedef tVector3<double> Vector3d;
#if defined(REAL_IS_FLOAT)
typedef Vector3f Vector3;
#else
typedef Vector3d Vector3;
#endif*/

template <class T>
inline void spherical_to_euclidian(T phi, T theta, tVector3<T>& d)
{
    d[0] = (T)cos(phi)*cos(theta);
    d[1] = (T)sin(phi)*cos(theta);
    d[2] = (T)sin(theta);
}


template <class T>
inline void euclidian_to_spherical(const tVector3<T>& d, T& phi, T& theta)
{
    theta = (T)(asin(d[2]));
    phi = (T)(atan2(d[1], d[0]));
}


template <class T>
inline T abs(T x)
{
	return (x >= 0) ? x : -x;
}


template <class T>
inline T clamp(T x, T low, T high)
{
	if (x < low) {
		return low;
	} else if (x > high) {
		return high;
	} else {
		return x;
	}
}
GTB_END_NAMESPACE


#ifndef OUTLINE
#include <gtb/math/math.ipp>
#endif

#include <gtb/math/amat.h>
#include <gtb/math/amatc.h>

#endif // GTB_MATH_INCLUDED
