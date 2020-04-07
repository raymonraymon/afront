
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


#ifndef GTB_LINE2_INCLUDED
#define GTB_LINE2_INCLUDED

#include <gtb/common.hpp>
#include <gtb/graphics/point2.hpp>
#include <gtb/graphics/vector2.hpp>

GTB_BEGIN_NAMESPACE

template <class T>
class tLine2 {
public:
    typedef T value_type;
    typedef tPoint2<T> Point2;

    tLine2();
    tLine2(const tLine2 &l);
    tLine2(const Point2 &p, const Point2 &q);
    tLine2(const Point2 &p, const tVector2<T> &d);

    tLine2 &operator=(const tLine2 &r);

    bool operator==(const tLine2 &l) const;
    bool operator!=(const tLine2 &l) const;

    tLine2 operator-() const;

    tVector2<T> &direction();
    const tVector2<T> &direction() const;
    Point2 &origin();
    const Point2 &origin() const;

    bool is_degenerate() const;


    // Returns an arbitrary point on the line.
    // point(i) == point(j), iff i == j.
    Point2 point(int i) const;

    // Returns a point on the line at time t
    Point2 point(value_type t) const;

    bool contains(const Point2 &p) const;
    Point2 projection(const Point2 &p) const;
    // distance between a point and the infinite line
    value_type distance(const Point2& p) const;
    tVector2<T> normal() const;

protected:
    Point2 _p;
    tVector2<T> _d;
};

GTB_GENERATE_CLASS_TYPEDEFS(Line2)

GTB_END_NAMESPACE

#ifndef OUTLINE
#include <gtb/graphics/line2.ipp>
#endif


#endif
