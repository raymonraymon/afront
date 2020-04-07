
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


#ifndef GTB_CONE_INCLUDED
#define GTB_CONE_INCLUDED

#include <gtb/common.hpp>
#include <gtb/graphics/point3.hpp>
#include <gtb/graphics/circle3.hpp>
#include <gtb/graphics/box3.hpp>


GTB_BEGIN_NAMESPACE


class Cone {
public:
	Cone();

	Cone(const Cone &cone);

	Cone(const Point3 &apex,
	     const Circle3 &base);

	Cone(const Point3 &apex,
	     const Point3 &base_center,
	     real_t base_radius);

	bool operator==(const Cone &c) const;
	bool operator!=(const Cone &c) const;

	const Point3 &apex() const;
	const Circle3 &base() const;
	real_t height() const;

	Box3 bounding_box() const;
	Point3 centroid() const;

	void render(unsigned num_sides = 6,
		    bool show_base = true) const;

protected:

	void compute_base_points(std::vector<Point3> &points,
				 unsigned num_sides) const;

	Point3 _apex;
	Circle3 _base;
};


GTB_END_NAMESPACE

#ifndef OUTLINE
#include <gtb/graphics/cone.ipp>
#endif

#endif // GTB_CONE_INCLUDED
