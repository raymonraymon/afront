
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
#include <gtb/graphics/view.hpp>
#include <gtb/graphics/line3.hpp>
#include <gtb/math/math.hpp>
#endif // WIN32


#ifdef OUTLINE
#define inline
#include <gtb/graphics/view.ipp>
#undef inline
#endif


using namespace std;


GTB_BEGIN_NAMESPACE


template <class T>
tView<T>::tView()
{
}

template <class T>
tView<T>::tView(const Camera &cam,
	   const Viewport &vp)
	: _camera(cam),
	  _viewport(vp)
{
}

template <class T>
tView<T>::tView(const tView &view)
	: _camera(view._camera),
	  _viewport(view._viewport)
{
}


template <class T>
void tView<T>::rotate_world(const Point3 &origin, value_type dx, value_type dy)
{
	if ((dx == 0) && (dy == 0)) {
		return;
	}
	// find angle
	value_type vx = (value_type) dx / (value_type) _viewport.width();
	value_type vy = (value_type) -dy / (value_type) _viewport.height();
	value_type theta = 2.0 * M_PI * (fabs(vx) + fabs(vy));

	// find axis
	Vector3 v = vx * _camera.right() + vy * _camera.up();
	Vector3 axis = _camera.towards().cross(v);
	tLine3<T> l(origin, axis);

	_camera.rotate(l, theta);
}


template <class T>
void tView<T>::translate_world(const Point3 &origin, value_type dx, value_type dy)
{
	if ((dx == 0) && (dy == 0)) {
		return;
	}
	value_type l = Point3::distance(origin, _camera.origin()) *
		tan(_camera.x_fov());
	value_type tx = -l * dx / _viewport.width();
	value_type ty = l * dy / _viewport.height();
	_camera.translate((tx * _camera.right()) + (ty * _camera.up()));
}


template <class T>
void tView<T>::scale_world(const Point3 &origin, value_type dx, value_type dy)
{
	(void) dx;
	if (dy == 0) {
		return;
	}
	value_type t = 2.0 * (value_type) dy / (value_type) _viewport.width();
	_camera.translate(t * (origin - _camera.origin()));
}


template <class T>
void tView<T>::scale_world(value_type factor)
{
	_camera.translate(factor * _camera.towards());
}


template <class T>
void tView<T>::load_camera(const vector<Light> &headlights) const
{
	_camera.load(_viewport, headlights);
}


template <class T>
void tView<T>::load_camera() const
{
	_camera.load(_viewport);
}

template <class T>
void tView<T>::gl_load_camera() const
{
	_camera.load(_viewport);
}


template <class T>
typename tView<T>::Point3 tView<T>::viewport_point(int x, int y) const
{
	Point2 vc(_viewport.center());
	value_type dx = (value_type) (x - vc.x()) / (value_type) (_viewport.width() / 2);
	value_type dy = (value_type) (y - vc.y()) / (value_type) (_viewport.height() / 2);
	Point3 wc = (_camera.origin() +
		     _camera.near_distance() * _camera.towards());
	Vector3 r =  _camera.right();
        r *= (_camera.near_distance() * tan(_camera.x_fov()));
	Vector3 u = _camera.up();
	u *= (_camera.near_distance() * tan(_camera.y_fov()));
	return wc + (dx * r) + (dy * u);
}


template <class T>
typename tView<T>::Ray3 tView<T>::world_ray(int x, int y) const
{
    return Ray3(_camera.origin(), viewport_point(x, y));
}


template <class T>
typename tView<T>::Point2 tView<T>::viewport_point(const Point3 &p) const
{
    // equivalent to gluProject
    Point3 q = _camera.world_to_camera(p);
    assert(!real::is_zero(q.z()));
    value_type dx = q.x() / (-q.z() * tan(_camera.x_fov()));
    value_type dy = q.y() / (-q.z() * tan(_camera.y_fov()));
    value_type x = _viewport.x_center() + dx * _viewport.width() * 0.5;
    value_type y = _viewport.y_center() + dy * _viewport.height() * 0.5;
    return Point2(x, y);
}

template <class T>
void tView<T>::set_viewport(const Viewport &vp)
{
    _viewport = vp;
    value_type fovx = _camera.x_fov();
    value_type w = vp.width();
    value_type h = vp.height();
    value_type fovy = atan((double)h/(double)w*tan(fovx));
    _camera.reset_fov(fovx, fovy);
}

template class tView<float>;
template class tView<double>;

GTB_END_NAMESPACE
