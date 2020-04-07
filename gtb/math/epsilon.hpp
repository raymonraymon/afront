
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



#ifndef GTB_EPSILON_INCLUDED
#define GTB_EPSILON_INCLUDED

#ifndef GTB_NAMESPACE_INCLUDED
#include <gtb/common.hpp>
#endif

GTB_BEGIN_NAMESPACE


bool eps_is_zero(double x, double eps);
bool eps_is_positive(double x, double eps);
bool eps_is_positive_or_zero(double x, double eps);
bool eps_is_negative(double x, double eps);
bool eps_is_negative_or_zero(double x, double eps);

bool eps_is_equal(double x, double y, double eps);
bool eps_is_greater(double x, double y, double eps);
bool eps_is_greater_or_equal(double x, double y, double eps);
bool eps_is_less(double x, double y, double eps);
bool eps_is_less_or_equal(double x, double y, double eps);

template<class T> T fp_precision();
template<> inline float fp_precision<float>() { return 1e-8f; }
template<> inline double fp_precision<double>() { return 1e-16; }

template<class T> T fp_sqrt_precision();
template<> inline float fp_sqrt_precision<float>() { return 1e-4f; }
template<> inline double fp_sqrt_precision<double>() { return 1e-16f; }

inline bool eps_is_zero(double x, double eps)
{
	return (x > -eps) && (x < eps);
}


inline bool eps_is_positive(double x, double eps)
{
	return x >= eps;
}


inline bool eps_is_positive_or_zero(double x, double eps)
{
	return x > -eps;
}


inline bool eps_is_negative(double x, double eps)
{
	return (x <= -eps);
}


inline bool eps_is_negative_or_zero(double x, double eps)
{
	return x < eps;
}


inline bool eps_is_equal(double x, double y, double eps)
{
	return eps_is_zero(x - y, eps);
}


inline bool eps_is_greater(double x, double y, double eps)
{
	return eps_is_positive(x - y, eps);
}


inline bool eps_is_greater_or_equal(double x, double y, double eps)
{
	return eps_is_positive_or_zero(x - y, eps);
}


inline bool eps_is_less(double x, double y, double eps)
{
	return eps_is_negative(x - y, eps);
}


inline bool eps_is_less_or_equal(double x, double y, double eps)
{
	return eps_is_negative_or_zero(x - y, eps);
}


GTB_END_NAMESPACE

#endif // GTB_EPSILON_INCLUDED
