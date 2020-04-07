
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
#include <gtb/gtb.hpp>
#include <gtb/graphics/cone.hpp>
#include <gtb/graphics/ogltools.h>
#endif // WIN32


#ifdef OUTLINE
#define inline
#include <gtb/graphics/cone.ipp>
#undef inline
#endif


using namespace std;


GTB_BEGIN_NAMESPACE


Cone::Cone()
	: _apex(0.0, 1.0, 0.0),
	  _base(Point3(0.0, -1.0, 0.0),
		Vector3(0.0, -1.0, 0.0),
		1.0)
{
}


Cone::Cone(const Cone &cone)
	: _apex(cone._apex),
	  _base(cone._base)
{
}


Cone::Cone(const Point3 &arg_apex,
	   const Circle3 &arg_base)
	: _apex(arg_apex),
	  _base(arg_base)
{
}


Cone::Cone(const Point3 &arg_apex,
	   const Point3 &arg_base_center,
	   real_t arg_base_radius)
	: _apex(arg_apex),
	  _base(arg_base_center,
		arg_base_center - arg_apex,
		arg_base_radius)
{
}


void Cone::render(unsigned num_sides,
		  bool show_base) const
{
	assert(num_sides >= 3);

	vector<Point3> points(num_sides);
	compute_base_points(points, num_sides);

	if (show_base) {
		_base.normal().load_as_normal();
		glBegin(GL_TRIANGLE_FAN);
		_base.center().load();
		for (unsigned i = 0; i < num_sides; i++) {
			points[i].load();
		}
		points[0].load();
		glEnd();
	}

	glBegin(GL_TRIANGLES);
	for (unsigned i = 0; i < num_sides; i++) {
		Vector3 n = Point3::normal(points[i],
					   _apex,
					   points[(i + 1) % num_sides]);
		n.load_as_normal();
		points[i].load();
		_apex.load();
		points[(i + 1) % num_sides].load();
	}
	glEnd();
}


Box3 Cone::bounding_box() const
{
	vector<Point3> points(5);
	compute_base_points(points, 4);
	points[4] = _apex;
	return Box3::bounding_box(points);
}


Point3 Cone::centroid() const
{
	return bounding_box().centroid();
}


void Cone::compute_base_points(vector<Point3> &points,
			       unsigned num_sides) const
{
	assert(points.size() >= num_sides);
    Matrix4 r = Vector3::rotation(Vector3::VECTOR3_POSITIVE_Z, _base.normal());
	Vector3 t(_base.center().x(), _base.center().y(), _base.center().z());
	for (unsigned i = 0; i < num_sides; i++) {
		real_t theta = ((real_t) i / (real_t) num_sides) * 2.0 * M_PI;
		points[i].reset(_base.radius() * cos(theta),
				_base.radius() * sin(theta),
				0.0);
		points[i] *= r;
		points[i] += t;
	}
}


GTB_END_NAMESPACE
