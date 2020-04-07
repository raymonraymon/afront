
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


#ifndef __MATHTOOLS_H
#define __MATHTOOLS_H

#include <gtb/common.hpp>
//#include <stdlib.h>

/*
 * Routines:
 *
 *  ipow - compute the integer power for small powers like x^2, x^3
 *  sqr - is just x*x
 *  max2, min2: same as std::min,max
 *  max3, min3: compute the maximal value of 3 values (optionally returning their index)
 *  assign_if_less -  if (x < y) x = y;
 *  logn - log2 of an integer - the number of bits required to represent the integer
 *  pwr2(x) returns 2^x
 *  bilinear - interporation of 4 numbers
 *  swap3 - same as std::swap, but cyclic for 3 values
 */
GTB_BEGIN_NAMESPACE

// Need to be place in the right header file

unsigned factorial(unsigned x);


template <class T> T ipow(T x, int exponent)
{
    T r=1;
    while (exponent-- > 0)
        r *= x;

    return r;
}

template <class T> T sqr(const T& x)
{
    return x*x;
}

template <class T> T max2(const T& x, const T& y)
{
    if (x > y)
    {
        return x;
    } else return y;
}

template <class T> T min2(const T& x, const T& y)
{
    if (x < y)
    {
        return x;
    } else return y;
}

template <class T> T max3(const T& x, const T& y, const T& z)
{
    if (x > y)
    {
        if (x > z) return x;
        else return z;
    } else if (y > z) return y;
    else return z;
}

/*
 * Same as max3, return the index [0..2] of the largest item.
 */
template <class T> T max3(const T& x, const T& y, const T& z, int& k)
{
    if (x > y)
    {
        if (x > z) {
			k = 0;
			return x;
		}
        else 
		{
			k = 2;
			return z;
		}
    }
	else if (y > z)
	{
		k = 1;
		return y;
	}
    else
	{
		k = 2;
		return z;
	}
}

template <class T> T min3(const T& x, const T& y, const T& z)
{
    if (x < y)
    {
        if (x < z) return x;
        else return z;
    } else if (y < z) return y;
    else return z;
}

/*
 * Same as max3, return the index [0..2] of the largest item.
 */
template <class T> T min3(const T& x, const T& y, const T& z, int& k)
{
    if (x < y)
    {
        if (x < z) {
			k = 0;
			return x;
		}
        else 
		{
			k = 2;
			return z;
		}
    }
	else if (y < z)
	{
		k = 1;
		return y;
	}
    else
	{
		k = 2;
		return z;
	}
}

template<class T>
void assign_if_less(T& x, const T& y)
{
    if (x < y) x = y;
}


int rand32();

/*
 * size_t logn(T v) - returns the number of bits required to store v
 */
template<class T>
inline unsigned logn(T v)
{
    T result = 0;

    while (v)
    {
        v >>= 1;
        ++result;
    }
    return result;
}

/*
 * T pwr2(T v) - return v to the power of 2
 */
template<class T>
inline T pwr2(T v)
{
    T result = 1;

    while (v)
    {
        result <<= 1;
        --v;
    }
    return result;
}

/*
 * T round2up(T v) - returns v rounded up to the nearst power of 2template<class T>
 */
template<class T>
inline T round2up(T v)
{
    return pwr2(logn(v));
}


/*
 * Bilinear interpolation
 *
 * Input:
 *
 *  v3(0,1)     v4(1,1)
 *       *----------*
 *       |          |
 *       |          |
 *       |          |
 *       |          |
 *       |  *(x,y)  |
 *       |          |
 *       *----------*
 *    v1(0,0)   v2 (1,0)
 *
 *
 * output: bilinear interpolation of v1..v4 to x,y.
 *
 * Assume: x,y are in the range of [0,1]
 *
 * Template argument REAL is float, double...
 */
template<class VAL, class REAL>
VAL bilinear(const VAL& v1, const VAL& v2, const VAL& v3, const VAL& v4, REAL x, REAL y)
{
    VAL v1_2_interp = (v1*(1-x) + v2*x);
    VAL v3_4_interp = (v3*(1-x) + v4*x);
	return ( v1_2_interp*(1-y) + v3_4_interp*y );
}

/*
 * Swap three values:
 *   v1 = v2
 *   v2 = v3
 *   v3 = v1
 */
template<class T>
void swap3(T& v1, T& v2, T& v3)
{
    T t = v1;
    v1 = v2;
    v2 = v3;
    v3 = v1;
}

//
// My abs function
//
template <class T>
inline T absv(T val)
{
    if (val < 0) return -val;
    else return val;
}

// nth falling power
inline double falling(double a, int n)
{
    double result = 1.0;
    for (int i = 0; i<n; ++i)
	result *= (double)(a-i);
    return result;
}

inline double ffactorial(int b)
{
    double result = 1.0;
    for (int i=1; i<=b; ++i)
	result *= (double)i;
    return result;
}

inline double fchoose(double a, int b)
{
    return falling(a,b) / ffactorial(b);
}

inline int choose(int a, int b)
{
    int result = 1;
    for (int i=0; i<b; i++) {
	result *= a-i;
	result /= i+1;
    }
    return result;
}

template<> inline double absv(double val) { return (val<0 ? -val : val); /*fabs(val)*/; }
template<> inline float absv(float val) { return (val<0 ? -val : val); /*fabsf(val)*/; }

GTB_END_NAMESPACE

#ifndef OUTLINE
#include <gtb/math/mathtools.inl>
#endif

#endif //__MATHTOOLS_H
