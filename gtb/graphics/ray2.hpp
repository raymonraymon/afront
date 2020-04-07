
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


#ifndef GTB_RAY2_INCLUDED
#define GTB_RAY2_INCLUDED

#include <gtb/common.hpp>
#include <gtb/real/real.hpp>
#include <gtb/graphics/point2.hpp>
#include <gtb/graphics/vector2.hpp>
#include <gtb/graphics/segment2.hpp>
#include <list>


GTB_BEGIN_NAMESPACE


//class Line2;
//class Polygon2;
//class Box2;

template<class T>
class tRay2 {
public:
    typedef T value_type;
    tRay2();
    tRay2(const tRay2 &r);
    tRay2(const tPoint2<T> &p, const tVector2<T> &d);
    tRay2(const tPoint2<T> &p, const tPoint2<T> &q);
    tRay2 &operator=(const tRay2 &r);

    bool operator==(const tRay2 &r) const;
    bool operator!=(const tRay2 &r) const;

    tRay2 operator-() const;

    const tPoint2<T> &source() const;
    const tVector2<T> &direction() const;
    bool is_degenerate() const;

    const tPoint2<T> point(value_type t) const;

    // Parametric value of closest point on ray.
    value_type t(const tPoint2<T> &p) const;

    tLine2<T> supporting_line() const;
    bool contains(const tPoint2<T> &p) const;
    tPoint2<T> projection(const tPoint2<T> &p) const;

    bool intersects(const tSegment2<T> &poly) const;

    bool intersects(const tSegment2<T> &poly,
                    value_type &hit_time) const;

    bool intersects(const tSegment2<T> &poly,
                    value_type &hit_time,
                    tPoint2<T> &hit_point) const;

    bool intersects(const tBox2<T> &box) const;

    bool intersects(const tBox2<T> &box, value_type &t1, value_type &t2) const;

    static const tRay2 POSITIVE_X;
    static const tRay2 NEGATIVE_X;
    static const tRay2 POSITIVE_Y;
    static const tRay2 NEGATIVE_Y;

protected:
    tPoint2<T> _source;
    tVector2<T> _direction;
};

GTB_GENERATE_CLASS_TYPEDEFS(Ray2)

GTB_END_NAMESPACE

#ifndef OUTLINE
#include <gtb/graphics/ray2.ipp>
#endif

#endif // GTB_RAY2_INCLUDED
