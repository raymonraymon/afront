
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
#include <gtb/graphics/ogltools.h>
#include <gtb/graphics/polygon3.hpp>
#include <gtb/graphics/segment3.hpp>
#endif // WIN32


#ifdef OUTLINE
#define inline
#include <gtb/graphics/polygon3.ipp>
#undef inline
#endif


using namespace std;


GTB_BEGIN_NAMESPACE


template<class T>
tPolygon3<T>::tPolygon3()
{
}


template<class T>
tPolygon3<T>::tPolygon3(const tPolygon3 &p)
{
	for (unsigned i = 0; i < p._points.size(); i++) {
		_points.push_back(p._points[i]);
	}
	_normal = p._normal;
	assert(_points.size() == p._points.size());
}


template<class T>
tPolygon3<T>::tPolygon3(const tPoint3<T> &a, const tPoint3<T> &b, const tPoint3<T> &c)
{
	_points.push_back(a);
	_points.push_back(b);
	_points.push_back(c);
	_normal = tPoint3<T>::normal(a, b, c);
}


template<class T>
tPolygon3<T>::tPolygon3(const tPoint3<T> &a,
		   const tPoint3<T> &b,
		   const tPoint3<T> &c,
		   const tVector3<T> &n)
{
	_points.push_back(a);
	_points.push_back(b);
	_points.push_back(c);
	_normal = n;
}


template<class T>
tPolygon3<T>::tPolygon3(const tPoint3<T> &a,
		   const tPoint3<T> &b,
		   const tPoint3<T> &c,
		   const tPoint3<T> &d)
{
	_points.push_back(a);
	_points.push_back(b);
	_points.push_back(c);
	_points.push_back(d);
	_normal = tPoint3<T>::normal(a, b, c);
}


template<class T>
tPolygon3<T>::tPolygon3(const tPoint3<T> &a,
		   const tPoint3<T> &b,
		   const tPoint3<T> &c,
		   const tPoint3<T> &d,
		   const tVector3<T> &n)
{
	_points.push_back(a);
	_points.push_back(b);
	_points.push_back(c);
	_points.push_back(d);
	_normal = n;
}


template<class T>
tPolygon3<T>::tPolygon3(const vector<tPoint3<T> > &points)
{
	assert(points.size() >= 3);
	for (unsigned i = 0; i < points.size(); i++) {
		_points.push_back(points[i]);
	}
	assert(_points.size() == points.size());
	_normal = tPoint3<T>::normal(_points[0], _points[1], _points[2]);
}


template<class T>
tPolygon3<T> &tPolygon3<T>::operator=(const tPolygon3 &p)
{
	if (&p != this) {
		_points.clear();
		for (unsigned i = 0; i < p._points.size(); i++) {
			_points.push_back(p._points[i]);
		}
		assert(_points.size() == p._points.size());
		_normal = p._normal;
	}
	return *this;
}


template<class T>
bool tPolygon3<T>::contains(const tPoint3<T> &p) const
{
	const tVector3<T> &n = normal();
	for (unsigned i = 0; i < _points.size(); i++) {
		unsigned j = (i + 1) % _points.size();
		const tPoint3<T> &pi = _points[i];
		const tPoint3<T> &pj = _points[j];
		tVector3<T> c = (pj - pi).cross(p - pi);
		if (real::is_negative(c.dot(n))) {
			return false;
		}
	}
	return true;
}


template<class T>
void tPolygon3<T>::render() const
{
	assert(_points.size() >= 3);
	glBegin(GL_POLYGON);
	for (unsigned i = 0; i < _points.size(); i++) {
		_points[i].load();
	}
	glEnd();
}

template class tPolygon3<float>;
template class tPolygon3<double>;

GTB_END_NAMESPACE
