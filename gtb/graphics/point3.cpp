
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
#include <gtb/graphics/point3.hpp>
#include <gtb/graphics/vector3.hpp>
#include <gtb/graphics/line3.hpp>
#endif // WIN32


#ifdef OUTLINE
#define inline
#include <gtb/graphics/point3.ipp>
#undef inline
#endif


using namespace std;


GTB_BEGIN_NAMESPACE

template<>
tPoint3<float> tPoint3<float>::ZERO(0,0,0);
template<>
tPoint3<double> tPoint3<double>::ZERO(0,0,0);


// const tPoint3 POINT3_ZERO(0.0, 0.0, 0.0);

template<class T>
tPoint3<T> &tPoint3<T>::scale(const tPoint3 &origin, value_type s)
{
	tVector3<T> t(origin._p);
	*this -= t;
	_p[0] *= s;
	_p[1] *= s;
	_p[2] *= s;
	*this += t;
	return *this;
}

template<class T>
tPoint3<T> &tPoint3<T>::rotate(const tVector3<T> &axis, value_type theta)
{
	tVector3<T> v(_p);
	v.rotate(axis, theta);
	this->reset(v.x(), v.y(), v.z());
	return *this;
}

template<class T>
tPoint3<T> &tPoint3<T>::rotate(const tLine3<T> &l, value_type theta)
{
	tPoint3 q = l.point(0);
	tVector3<T> v = *this - q;
	v.rotate(l.direction(), theta);
	*this = q + v;
	return *this;
}

template<class T>
tPoint3<T> tPoint3<T>::centroid(const vector<tPoint3> &v)
{
	value_type cx = 0.0;
	value_type cy = 0.0;
	value_type cz = 0.0;

	for (unsigned i = 0; i < v.size(); i++) {
		cx += v[i].x() / v.size();
		cy += v[i].y() / v.size();
		cz += v[i].z() / v.size();
	}
	return tPoint3(cx, cy, cz);
}

template<>
void tPoint3<double>::read(FILE *fp)
{
	read_double(&_p[0], fp);
	read_double(&_p[1], fp);
	read_double(&_p[2], fp);
}

template<>
 void tPoint3<double>::write(FILE *fp) const
{
	write_double(_p[0], fp);
	write_double(_p[1], fp);
	write_double(_p[2], fp);
}

template<>
 void tPoint3<float>::read(FILE *fp)
{
	read_float(&_p[0], fp);
	read_float(&_p[1], fp);
	read_float(&_p[2], fp);
}

template<>
void tPoint3<float>::write(FILE *fp) const
{
	write_float(_p[0], fp);
	write_float(_p[1], fp);
	write_float(_p[2], fp);
}


template class tPoint3<float>;
template class tPoint3<double>;

GTB_END_NAMESPACE
