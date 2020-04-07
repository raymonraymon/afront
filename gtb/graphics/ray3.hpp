
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


#ifndef GTB_RAY3_INCLUDED
#define GTB_RAY3_INCLUDED

#include <gtb/common.hpp>
#include <gtb/real/real.hpp>
#include <gtb/graphics/point3.hpp>
#include <gtb/graphics/vector3.hpp>
#include <gtb/graphics/polygon3.hpp>
#include <list>


GTB_BEGIN_NAMESPACE


template <class T> class tLine3;
template <class T> class tPolygon3;
template <class T> class tBox3;

template<class T>
class tRay3 {
public:
    typedef T value_type;
    tRay3();
    tRay3(const tRay3 &r);
    tRay3(const tPoint3<T> &p, const tVector3<T> &d);
    tRay3(const tPoint3<T> &p, const tPoint3<T> &q);
    tRay3 &operator=(const tRay3 &r);

    bool operator==(const tRay3 &r) const;
    bool operator!=(const tRay3 &r) const;

    tRay3 operator-() const;

    const tPoint3<T> &source() const;
    const tVector3<T> &direction() const;
    bool is_degenerate() const;

    const tPoint3<T> point(value_type t) const;

    // Parametric value of closest point on ray.
    value_type t(const tPoint3<T> &p) const;

    tLine3<T> supporting_line() const;
    bool contains(const tPoint3<T> &p) const;
    tPoint3<T> projection(const tPoint3<T> &p) const;

    bool intersects(const tPolygon3<T> &poly) const;

    bool intersects(const tPolygon3<T> &poly,
                    value_type &hit_time) const;

    bool intersects(const tPolygon3<T> &poly,
                    value_type &hit_time,
                    tPoint3<T> &hit_point) const;

    bool intersects(const tBox3<T> &box) const;

    bool intersects(const tBox3<T> &box, value_type &t1, value_type &t2) const;

    static const tRay3 POSITIVE_X;
    static const tRay3 NEGATIVE_X;
    static const tRay3 POSITIVE_Y;
    static const tRay3 NEGATIVE_Y;
    static const tRay3 POSITIVE_Z;
    static const tRay3 NEGATIVE_Z;

protected:
    tPoint3<T> _source;
    tVector3<T> _direction;
};

GTB_GENERATE_CLASS_TYPEDEFS(Ray3)

GTB_END_NAMESPACE

#ifndef OUTLINE
#include <gtb/graphics/ray3.ipp>
#endif

#endif // GTB_RAY3_INCLUDED
