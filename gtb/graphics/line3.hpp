
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


#ifndef GTB_LINE3_INCLUDED
#define GTB_LINE3_INCLUDED

#include <gtb/common.hpp>
#include <gtb/graphics/point3.hpp>
#include <gtb/graphics/vector3.hpp>

GTB_BEGIN_NAMESPACE

#if 0
template<class T>
class tRay3;
#endif

template<class T>
class tLine3 {
public:
    typedef T value_type;

    typedef tPoint3<T> Point3;

    tLine3();
    tLine3(const tLine3 &l);
    tLine3(const Point3 &p, const Point3 &q);
    tLine3(const Point3 &p, const tVector3<T> &d);
//    explicit tLine3(const tRay3<T> &r);
    tLine3 &operator=(const tLine3 &r);
    
    bool operator==(const tLine3 &l) const;
    bool operator!=(const tLine3 &l) const;
    
    tLine3 operator-() const;
    
    tVector3<T>& direction();
    const tVector3<T> &direction() const;
    Point3& origin();
    const Point3& origin() const;

    bool is_degenerate() const;

    // Returns an arbitrary point on the line.
    // point(i) == point(j), iff i == j.
    Point3 point(int i) const;

    //
    // Return a point on the line at time t
    //
    Point3 point(value_type t) const;

    bool contains(const Point3 &p) const;
    Point3 projection(const Point3 &p) const;
    // distance between a point and the infinite line
    value_type distance(const Point3& p) const;

protected:
    Point3 _p;
    tVector3<T> _d;

#if 0
    friend bool LineDistance(const tLine3<T>& L1, const tLine3<T>& L2, Point3& pc, Point3& qc);
    friend bool LineDistance(const tLine3<T>& L1, const tLine3<T>& L2, T& sc, T& tc);
#endif
};

#if 0
extern const tLine3 LINE3_POSITIVE_X;
extern const tLine3 LINE3_NEGATIVE_X;
extern const tLine3 LINE3_POSITIVE_Y;
extern const tLine3 LINE3_NEGATIVE_Y;
extern const tLine3 LINE3_POSITIVE_Z;
extern const tLine3 LINE3_NEGATIVE_Z;
#endif

/*
 * Return two points pc on L1 and qc on L2 that are the closest
 * points from L1 to L2.
 * ||pc-qc|| is the distnace between L1 and L2
 * 
 * If the lines are parallel, return false, otherwise return true.
 *
 * Assume the line directions are normalized
 */
template<class T>
bool LineDistance(const tLine3<T>& L1, const tLine3<T>& L2, tPoint3<T>& pc, tPoint3<T>& qc);

/*
 * Return two points pc on L1 and qc on L2 that are the closest
 * points from L1 to L2.
 * "sc" and "tc" are the time of the closest points in L1, L2 respectively
 * 
 * If the lines are parallel, return false, otherwise return true.
 *
 * Assume the line directions are normalized
 */
template<class T>
bool LineDistance(const tLine3<T>& L1, const tLine3<T>& L2, T& sc, T& tc);

typedef tLine3<float> Line3f;
typedef tLine3<double> Line3d;
#if defined(REAL_IS_FLOAT)
typedef Line3f Line3;
#else
typedef Line3d Line3;
#endif

GTB_END_NAMESPACE

#ifndef OUTLINE
#include <gtb/graphics/line3.ipp>
#endif

#endif // GTB_LINE3_INCLUDED
