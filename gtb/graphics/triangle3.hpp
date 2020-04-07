
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


#ifndef GTB_TRIANGLE3_INCLUDED
#define GTB_TRIANGLE3_INCLUDED

#include <gtb/common.hpp>
#include <gtb/graphics/point3.hpp>
#include <gtb/graphics/point2.hpp>
#include <gtb/graphics/box3.hpp>
#include <gtb/graphics/circle3.hpp>


GTB_BEGIN_NAMESPACE


//class Vector3;
//class Circle3;

template <class T>
class tTriangle3 {
public:
    typedef T value_type;
    typedef tPoint2<T> Point2;
    typedef tPoint3<T> Point3;
    typedef tVector3<T> Vector3;
    typedef tCircle3<T> Circle3;
    typedef tBox3<T> Box3;

    tTriangle3();
    tTriangle3(const tTriangle3 &t);
    tTriangle3(const Point3 &A, const Point3 &B, const Point3 &C);
    tTriangle3 &operator=(const tTriangle3 &t);

    bool operator==(const tTriangle3 &t) const;
    bool operator!=(const tTriangle3 &t) const;

    const Point3 &A() const;
    const Point3 &get_A() const;
    const Point3 &B() const;
    const Point3 &get_B() const;
    const Point3 &C() const;
    const Point3 &get_C() const;

	const Point3 &operator[](int i) const;


    Vector3 normal() const;
    value_type area() const;
    Circle3 circumcircle() const;
    Box3 bounding_box() const;
	Point3 closest_point(const Point3 &from) const;
    void subdivide(value_type max_area, std::vector<tTriangle3> &v) const;
    void get_barycentric_coordinates(const Point3 &p,
                                     value_type &wa,
                                     value_type &wb,
                                     value_type &wc) const;

	// t returns the distance along the line - may be negative
	bool intersect_ray(const Point3 &p0, const Vector3 &d, value_type &t) const;

	// t returns the fraction along the line - 0=p1, 1=p2
	bool intersect_segment(const Point3 &p1, const Point3 &p2, value_type &t) const;

	// get segment where two triangles intersect
	static bool Intersection(const tTriangle3<T> &tri1, const tTriangle3<T> &tri2, typename tTriangle3<T>::Point3 &i1, typename tTriangle3<T>::Point3 &i2);

protected:
    Point3 m_A, m_B, m_C;
};






typedef tTriangle3<float> Triangle3f;
typedef tTriangle3<double> Triangle3d;
#if defined(REAL_IS_FLOAT)
typedef Triangle3f Triangle3;
#else
typedef Triangle3d Triangle3;
#endif

GTB_END_NAMESPACE

#ifndef OUTLINE
#include <gtb/graphics/triangle3.ipp>
#endif

#endif // GTB_TRIANGLE3_INCLUDED
