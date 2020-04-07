
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
#include <gtb/graphics/ray2.hpp>
#include <gtb/graphics/segment2.hpp>
#include <gtb/graphics/box2.hpp>
#endif // WIN32


#ifdef OUTLINE
#define inline
#include <gtb/graphics/ray2.ipp>
#undef inline
#endif


using namespace std;


GTB_BEGIN_NAMESPACE


template<>
const tRay2<float> tRay2<float>::POSITIVE_X(tPoint2<float>(0.0f, 0.0f), tVector2<float>(1.0f, 0.0f ));
template<>
const tRay2<float> tRay2<float>::NEGATIVE_X(tPoint2<float>(0.0f, 0.0f), tVector2<float>(-1.0f, 0.0f));
template<>
const tRay2<float> tRay2<float>::POSITIVE_Y(tPoint2<float>(0.0f, 0.0f), tVector2<float>(0.0f, 1.0f ));
template<>
const tRay2<float> tRay2<float>::NEGATIVE_Y(tPoint2<float>(0.0f, 0.0f), tVector2<float>(0.0f, -1.0f));

template<>
const tRay2<double> tRay2<double>::POSITIVE_X(tPoint2<double>(0.0, 0.0), tVector2<double>(1.0, 0.0 ));
template<>
const tRay2<double> tRay2<double>::NEGATIVE_X(tPoint2<double>(0.0, 0.0), tVector2<double>(-1.0, 0.0));
template<>
const tRay2<double> tRay2<double>::POSITIVE_Y(tPoint2<double>(0.0, 0.0), tVector2<double>(0.0, 1.0 ));
template<>
const tRay2<double> tRay2<double>::NEGATIVE_Y(tPoint2<double>(0.0, 0.0), tVector2<double>(0.0, -1.0));


template<class T>
tRay2<T>::tRay2()
	: _source(tPoint2<T>((T)0.0, (T)0.0)),
	  _direction(tVector2<T>((T)1.0, (T)0.0))
{
}


template<class T>
tRay2<T>::tRay2(const tRay2 &r)
	: _source(r._source),
	  _direction(r._direction)
{
}


template<class T>
tRay2<T>::tRay2(const tPoint2<T> &p, const tVector2<T> &d)
	: _source(p),
	  _direction(d)
{
	_direction.normalize();
}


template<class T>
tRay2<T>::tRay2(const tPoint2<T> &p, const tPoint2<T> &q)
	: _source(p),
	  _direction(q - p)
{
	_direction.normalize();
}


template<class T>
bool tRay2<T>::intersects(const tSegment2<T> &poly) const
{
	value_type hit_time;
	tPoint2<T> hit_point;
	return intersects(poly, hit_time, hit_point);
}


template<class T>
bool tRay2<T>::intersects(const tSegment2<T> &poly, value_type &hit_time) const
{
	tPoint2<T> hit_point;
	return intersects(poly, hit_time, hit_point);
}


template<class T>
bool tRay2<T>::intersects(const tSegment2<T> &poly,
		      value_type &hit_time,
		      tPoint2<T> &hit_point) const
{
	const tVector2<T> &d = _direction;
	const tPoint2<T> &s = _source;
	value_type D = poly.D();
	tVector2<T> n = poly.normal();

	value_type nd = n * d;
	if (real::is_zero(nd)) {
		return false;
	}
	value_type ns = n * s;
	hit_time = -(D + ns) / nd;
	if (real::is_negative(hit_time)) {
		return false;
	}
	hit_point = point(hit_time);
	return poly.contains(hit_point);
}


template<class T>
bool tRay2<T>::intersects(const tBox2<T> &box) const
{
	value_type t1, t2;
	return intersects(box, t1, t2);
}


// See Moller and Haines, p. 301.
template<class T>
bool tRay2<T>::intersects(const tBox2<T> &box, value_type &t1, value_type &t2) const
{
	const tPoint2<T> &pmin = box.min_point();
	const tPoint2<T> &pmax = box.max_point();
	value_type tmin = -real::MAX;
	value_type tmax = real::MAX;
	tVector2<T> p = box.centroid() - _source;
	// for each slab (pair of parallel planes)
	for (unsigned i = 0; i < 2; i++) {
		value_type e = p[i];
		value_type f = _direction[i];
		value_type h = 0.5 * (pmax[i] - pmin[i]);
		if (!real::is_zero(f)) {
			// ray hits slab
			value_type finv = 1.0 / f;
			t1 = (e + h) * finv;
			t2 = (e - h) * finv;
			if (t1 > t2) {
				value_type tmp = t1;
				t1 = t2;
				t2 = tmp;
			}
			if (t1 > tmin) {
				tmin = t1;
			}
			if (t2 < tmax) {
				tmax = t2;
			}
			if (tmin > tmax) {
				// ray misses the box
				return false;
			}
			if (tmax < 0.0) {
				// box behind ray's origin
				return false;
			}
		} else {
			// ray parallel to slab
			if ((-e - h > 0.0) || (-e + h < 0.0)) {
				// ray outside slab
				return false;
			}
		}
	}
	assert(t2 >= 0.0);
	assert(t2 >= t1);
	// It's OK for t1 to be < 0 (ray origin inside the box).
	return true;
}

template class tRay2<float>;
template class tRay2<double>;

GTB_END_NAMESPACE
