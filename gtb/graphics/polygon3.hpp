
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


#ifndef GTB_POLYGON3_INCLUDED
#define GTB_POLYGON3_INCLUDED

#include <gtb/common.hpp>
#include <gtb/graphics/point3.hpp>
#include <gtb/graphics/vector3.hpp>
#include <gtb/graphics/box3.hpp>
#include <gtb/graphics/plane.hpp>

GTB_BEGIN_NAMESPACE


//class Plane;
//class tSegment3<T>;

template <class T> class tBox3;

template<class T>
class tPolygon3 {
public:
    typedef T value_type;
    tPolygon3();

    tPolygon3(const tPolygon3 &p);

    tPolygon3(const tPoint3<T> &a,
              const tPoint3<T> &b,
              const tPoint3<T> &c);

    tPolygon3(const tPoint3<T> &a,
              const tPoint3<T> &b,
              const tPoint3<T> &c,
              const tVector3<T> &normal);

    tPolygon3(const tPoint3<T> &a,
              const tPoint3<T> &b,
              const tPoint3<T> &c,
              const tPoint3<T> &d);

    tPolygon3(const tPoint3<T> &a,
              const tPoint3<T> &b,
              const tPoint3<T> &c,
              const tPoint3<T> &d,
              const tVector3<T> &normal);

    tPolygon3(const std::vector<tPoint3<T> > &points);

    tPolygon3 &operator=(const tPolygon3 &p);

    const tVector3<T> &normal() const;
    value_type D() const;
    unsigned num_points() const;
    const tPoint3<T> &point(unsigned i) const;
    tPlane<T> plane() const;

    tBox3<T> bounding_box() const;
    tPoint3<T> centroid() const;

    bool contains(const tPoint3<T> &p) const; // assumes convex, ccw

    bool intersects(const tSegment3<T> &segment) const; // assumes convex, ccw

    bool intersects(const tSegment3<T> &segment, // assumes convex, ccw
                    value_type &hit_time,
                    tPoint3<T> &hit_point) const;

    bool intersects(const tPolygon3 &p) const; // assumes convex, ccw

    void render() const;

protected:
    bool asymmetrically_intersects(const tPolygon3 &p) const;

    std::vector<tPoint3<T> > _points;
    tVector3<T> _normal;
};

typedef tPolygon3<float> Polygon3f;
typedef tPolygon3<double> Polygon3d;

#ifdef REAL_IS_FLOAT
typedef Polygon3f Polygon3;
#else
typedef Polygon3d Polygon3;
#endif


GTB_END_NAMESPACE

#ifndef OUTLINE
#include <gtb/graphics/polygon3.ipp>
#endif

#endif // GTB_POLYGON3_INCLUDED
