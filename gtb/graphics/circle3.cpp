
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
#include <gtb/graphics/circle3.hpp>
#include <gtb/math/math.hpp>
#include <gtb/graphics/ogltools.h>
#endif // WIN32


#ifdef OUTLINE
#define inline
#include <gtb/graphics/circle3.ipp>
#undef inline
#endif


GTB_BEGIN_NAMESPACE


template <class T>
tCircle3<T>::tCircle3()
    : _center(0.0, 0.0, 0.0),
      _normal(0.0, 0.0, 1.0),
      _radius(1.0)
{
}


template <class T>
tCircle3<T>::tCircle3(const Point3 &c,
                      const Vector3 &n,
                      value_type r)
    : _center(c),
      _normal(n),
      _radius(r)
{
    _normal.normalize();
}


template <class T>
tCircle3<T>::tCircle3(const tCircle3 &c)
    : _center(c._center),
      _normal(c._normal),
      _radius(c._radius)
{
    _normal.normalize();
}


template <class T>
tCircle3<T>::tCircle3(const Point3 &A,
                      const Point3 &B,
                      const Point3 &C)
{
    Vector3 AB = B - A;
    Vector3 BC = C - B;
    Vector3 CA = A - C;

    value_type a = BC.length();
    value_type b = CA.length();
    value_type c = AB.length();

    value_type cos_A = -CA.dot(AB) / (b * c);    
    value_type cos_B = -AB.dot(BC) / (c * a);
    value_type cos_C = -BC.dot(CA) / (a * b);

    value_type sin_A = sqrt(1.0 - cos_A * cos_A);
    value_type R = a / (2.0 * sin_A);
    value_type area = 0.5 * b * c * sin_A;

    value_type t1 = 0.5 * a * R * cos_A / area;
    value_type t2 = 0.5 * b * R * cos_B / area;
    value_type t3 = 0.5 * c * R * cos_C / area;

    _center.reset((t1 * A.x()) + (t2 * B.x()) + (t3 * C.x()),
                  (t1 * A.y()) + (t2 * B.y()) + (t3 * C.y()),
                  (t1 * A.z()) + (t2 * B.z()) + (t3 * C.z()));
    _normal = AB.cross(-CA);
    _normal.normalize();
    _radius = R;
}


template <class T>
void tCircle3<T>::render(unsigned num_sides) const
{
    tMatrix4<T> r = Vector3::rotation(Vector3::VECTOR3_POSITIVE_Z, _normal);
    Vector3 t(_center.x(), _center.y(), _center.z());
    _normal.load_as_normal();
    glBegin(GL_TRIANGLE_FAN);
    _center.load();
    for (unsigned i = 0; i <= num_sides; i++) {
        value_type theta = ((value_type) i / (value_type) num_sides) * 2.0 * M_PI;
        Point3 p(_radius * cos(theta),
                 _radius * sin(theta),
                 0.0);
        p *= r;
        p += t;
        p.load();
    }
    glEnd();
}


template class tCircle3<float>;
template class tCircle3<double>;


GTB_END_NAMESPACE
