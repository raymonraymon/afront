
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


#ifndef GTB_REAL_INCLUDED
#define GTB_REAL_INCLUDED

#include <gtb/common.hpp>
#include <float.h>
#include <stdio.h>
#include <gtb/io/io.hpp>

GTB_BEGIN_NAMESPACE


#if defined(REAL_IS_FLOAT)
typedef float real_t;
#elif defined(REAL_IS_DOUBLE)
typedef double real_t;
#else
#error "Bad type for real_t"
#endif


class real {
public:
	static real_t EPS;
	static const real_t MIN;
	static const real_t MAX;

	static bool is_zero(real_t x, real_t eps = EPS);
	static bool is_positive(real_t x, real_t eps = EPS);
	static bool is_positive_or_zero(real_t x, real_t eps = EPS);
	static bool is_negative(real_t x, real_t eps = EPS);
	static bool is_negative_or_zero(real_t x, real_t eps = EPS);

	static bool is_equal(real_t x, real_t y, real_t eps = EPS);
	static bool is_greater(real_t x, real_t y, real_t eps = EPS);
	static bool is_greater_or_equal(real_t x, real_t y, real_t eps = EPS);
	static bool is_less(real_t x, real_t y, real_t eps = EPS);
	static bool is_less_or_equal(real_t x, real_t y, real_t eps = EPS);

	static void set_eps(real_t eps);
	static float epsf();
	static double epsd();

	static int scan(FILE *fp, real_t *x);
};

template<class T>
class treal {
public:
    typedef T value_type;
    static value_type EPS;
//    static const value_type MIN;
//    static const value_type MAX;

    static bool is_zero(value_type x, value_type eps = EPS);
    static bool is_equal(value_type x, value_type y, value_type eps = EPS);
    static int scan(FILE *fp, value_type *x);
    static bool is_positive(value_type x);
};

template <class T> struct type_traits
{
    static void read(T *v, FILE *fp);
    static void write(T v, FILE *fp);
    static T cos(T v);
    static T sin(T v);
};

template <>
struct type_traits<float>
{
    static void read(float *f, FILE *fp) { read_float(f, fp); }
    static void write(float f, FILE *fp) { write_float(f, fp); }
//	static float cos(float v) { return cosf(v); }
//	static float sin(float v) { return sinf(v); }
};

template <>
struct type_traits<double>
{
    static void read(double *f, FILE *fp) { read_double(f, fp); }
    static void write(double f, FILE *fp) { write_double(f, fp); }
//	static double cos(double v) { return ::cos(v); }
//	static double sin(double v) { return ::sin(v); }
};

GTB_END_NAMESPACE


#ifndef OUTLINE
#include <gtb/real/real.ipp>
#endif


#endif // GTB_REAL_INCLUDED
