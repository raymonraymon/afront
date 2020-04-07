
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


#ifndef GTB_CIRCLE3_INCLUDED
#define GTB_CIRCLE3_INCLUDED

#include <gtb/common.hpp>
#include <gtb/graphics/point3.hpp>
#include <gtb/graphics/vector3.hpp>


GTB_BEGIN_NAMESPACE


template <class T>
class tCircle3 {
public:
    typedef T value_type;
    typedef tPoint3<T> Point3;
    typedef tVector3<T> Vector3;

    tCircle3();
    tCircle3(const tCircle3 &c);
    tCircle3(const Point3 &center,
             const Vector3 &normal,
             value_type radius);
    tCircle3(const Point3 &A, const Point3 &B, const Point3 &C);
    tCircle3 &operator=(const tCircle3 &c);

    bool operator==(const tCircle3 &c) const;
    bool operator!=(const tCircle3 &c) const;

    const Point3 &center() const;
    const Vector3 &normal() const;
    value_type radius() const;

    void render(unsigned num_sides = 6) const;

protected:
    Point3 _center;
    Vector3 _normal;
    value_type _radius;
};

typedef tCircle3<float> Circle3f;
typedef tCircle3<double> Circle3d;
#if defined(REAL_IS_FLOAT)
typedef Circle3f Circle3;
#else
typedef Circle3d Circle3;
#endif


GTB_END_NAMESPACE

#ifndef OUTLINE
#include <gtb/graphics/circle3.ipp>
#endif

#endif // GTB_CIRCLE3_INCLUDED
