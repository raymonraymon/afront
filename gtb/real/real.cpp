
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


#include <gtb/gtb.hpp>
#ifndef WIN32
#include <gtb/real/real.hpp>
#endif // WIN32


#ifdef OUTLINE
#define inline
#include <gtb/real/real.ipp>
#undef inline
#endif


GTB_BEGIN_NAMESPACE


#if defined(REAL_IS_FLOAT)
real_t real::EPS = real::epsf();
const real_t real::MIN = FLT_MIN;
const real_t real::MAX = FLT_MAX;
#elif defined(REAL_IS_DOUBLE)
real_t real::EPS = real::epsd();
const real_t real::MIN = DBL_MIN;
const real_t real::MAX = DBL_MAX;
#else
#error "Bad type for real_t"
#endif


template <class T>
static T epst()
{
	T a = 1.0;
	T e = 1.0;
	T f = e / 2.0;
	T b = a + f;

	while (a != b) {
		e = f;
		f /= 2.0;
		b = a + f;
	}
	return e;
}


float real::epsf()
{
	return epst<float>();
}


double real::epsd()
{
	return epst<double>();
}


int real::scan(FILE *fp, real_t *x)
{
#ifdef REAL_IS_FLOAT
	return fscanf(fp, "%f", x);
#else
	return fscanf(fp, "%lf", x);
#endif
}


template<>
treal<float>::value_type treal<float>::EPS = 1e-12f;
template<>
treal<double>::value_type treal<double>::EPS = 1e-24;


template<>
int treal<float>::scan(FILE *fp, value_type *x)
{
	return fscanf(fp, "%f", x);
}

template<>
int treal<double>::scan(FILE *fp, value_type *x)
{
	return fscanf(fp, "%lf", x);
}


template class treal<float>;
template class treal<double>;

GTB_END_NAMESPACE
