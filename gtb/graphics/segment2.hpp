
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


#ifndef GTB_SEGMENT2_INCLUDED
#define GTB_SEGMENT2_INCLUDED

#include <gtb/common.hpp>
#include <gtb/graphics/point2.hpp>
#include <gtb/graphics/line2.hpp>

GTB_BEGIN_NAMESPACE

template<class T> class tBox2;

//class Vector2;
//class Box2;
//class Line2;

template<class T>
class tSegment2 {
public:
    typedef T value_type;
    typedef tPoint2<T> Point2;
    typedef tVector2<T> Vector2;
    typedef tLine2<T> Line2;

    tSegment2();
    tSegment2(const tSegment2 &s);
    tSegment2(const Point2 &source, const Point2 &target);
    tSegment2 &operator=(const tSegment2 &s);

    bool operator==(const tSegment2 &s) const;
    bool operator!=(const tSegment2 &s) const;

    const Point2 &source() const;
    const Point2 &target() const;

    const Point2 &min() const;
    const Point2 &max() const;

    value_type D() const;
    value_type squared_length() const;
    value_type length() const;
    Vector2 direction() const;
    tLine2<T> supporting_line() const;
    tBox2<T> bounding_box() const;
    tVector2<T> normal() const;

    bool is_degenerate() const;
    bool contains(const Point2 &p) const;
    bool collinear_contains(const Point2 &p) const;

    tSegment2 operator-() const;

    // return the distance between the point and the line segment
    T distance(const Point2& p) const;
    T squared_distance(const Point2& p) const;

protected:
    Point2 _p, _q;
};

GTB_GENERATE_CLASS_TYPEDEFS(Segment2)


GTB_END_NAMESPACE

#ifndef OUTLINE
#include <gtb/graphics/segment2.ipp>
#endif

#endif // GTB_SEGMENT3_INCLUDED
