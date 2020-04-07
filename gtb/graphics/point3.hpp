
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


#ifndef GTB_POINT3_INCLUDED
#define GTB_POINT3_INCLUDED

#include <gtb/common.hpp>
#include <gtb/real/real.hpp>
//#include <gtb/math/math.hpp>
#include <vector>

GTB_BEGIN_NAMESPACE


template<class T>
class tVector3;

template<class T>
class tLine3;

template<class T>
class tMatrix4;

template<class T>
class tPoint3 {
public:
    typedef T value_type;
    tPoint3();
    template <class T2> explicit tPoint3(const tPoint3<T2> &p);
    tPoint3(const tVector3<T> &v);
    tPoint3(value_type px, value_type py, value_type pz);
    tPoint3(const value_type p[3]);
    tPoint3 &operator=(const tPoint3 &p);

    value_type x() const;
    value_type y() const;
    value_type z() const;
    value_type get_x() const;
    value_type get_y() const;
    value_type get_z() const;

    void set_x(value_type px);
    void set_y(value_type py);
    void set_z(value_type pz);
    tPoint3 &reset(value_type px, value_type py, value_type pz);

    value_type operator[](unsigned i) const;
    value_type &operator[](unsigned i);

    bool operator==(const tPoint3 &p) const;
    bool operator!=(const tPoint3 &p) const;
    bool is_zero() const;

    tPoint3 &operator*=(const tMatrix4<T> &m);
    tPoint3 &operator+=(const tVector3<T> &v);
    tPoint3 &operator-=(const tVector3<T> &v);

    void add(const tPoint3& p);
    void add_scaled(const tPoint3& p, value_type s);
    void scalar_scale(value_type s);
    tPoint3 scalar_scaled(value_type s) const;

    template<class T2>
    friend tPoint3<T2> operator*(const tMatrix4<T2> &m, const tPoint3<T2> &p);
    template<class T2>
    friend tPoint3<T2> operator+(const tPoint3<T2> &p, const tVector3<T2> &v);
    template<class T2>
    friend tVector3<T2> operator-(const tPoint3<T2> &p, const tPoint3<T2> &q);
    template<class T2>
    friend tPoint3<T2> operator-(const tPoint3<T2> &p, const tVector3<T2> &v);

    tPoint3 &scale(const tPoint3 &origin, value_type s);
    tPoint3 &translate(const tVector3<T> &v);
    tPoint3 &translate(value_type dx, value_type dy, value_type dz);
    tPoint3 &rotate(const tVector3<T> &axis, value_type theta);
    tPoint3 &rotate(const tLine3<T> &l, value_type theta);
    tPoint3 &transform(const tMatrix4<T> &m);
    tPoint3 &affine_transform(const tMatrix4<T> &m);

    static value_type distance(const tPoint3<T> &p, const tPoint3<T> &q);
    static value_type squared_distance(const tPoint3<T> &p, const tPoint3<T> &q);

    static bool collinear(const tPoint3 &A,
                          const tPoint3 &B,
                          const tPoint3 &C);

    static tVector3<T> normal(const tPoint3 &A,
                              const tPoint3 &B,
                              const tPoint3 &C);

    static tPoint3 midpoint(const tPoint3 &A, const tPoint3 &B);
    static tPoint3 inbetween(const tPoint3 &A, const tPoint3 &B, value_type alpha);

    static tPoint3 centroid(const tPoint3 &A,
                            const tPoint3 &B,
                            const tPoint3 &C);
#if 0
    static tPoint3 centroid(const tPoint3 &A,
                            const tPoint3 &B,
                            const tPoint3 &C,
                            const tPoint3 &D);
#endif
    static tPoint3 centroid(const std::vector<tPoint3> &v);

    void load() const;
    void gl_vertex() const;
    //void render() const;

    void read(FILE *fp);
    void write(FILE *fp) const;

    friend class tVector3<T>;

#if 0
    friend std::istream &operator>>(std::istream &s, tPoint3<T> &p);
    friend std::ostream &operator<<(std::ostream &s, const tPoint3<T> &p);
#endif

    static tPoint3 ZERO;

protected:
    value_type _p[3];
};

/*
 * Return > 0 if left-turn
 *        < 0 if right-turn
 *        = 0 if colinear
 * value = area of parallelogram that is defined by the 3 points
 */
template<class T>
inline T LeftTurn(const tPoint3<T>& p0, const tPoint3<T>& p1, const tPoint3<T>& p2) 
{
    return (p1.x()-p0.x())*(p2.y()-p1.y()) -
        (p2.x()-p1.x())*(p1.y()-p0.y());
}


/// by Olga - this function returns true if the given
/// point is inside a triangle defined by p0,p1,p2
/// Pre-condition: both point and *this are on xy plane
template<class T>
inline bool ContainsPoint(const tPoint3<T>& p0, const tPoint3<T>& p1, const tPoint3<T>& p2, const tPoint3<T>& point)
{
    T L0 = LeftTurn(p0, p1, point);
    T L1 = LeftTurn(p1, p2, point);
    T L2 = LeftTurn(p2, p0, point);

    if ((fabs(L0) < 1e-12) || (fabs(L1) < 1e-12) || (fabs(L2) < 1e-12)) return true;

    return ((L0<=0 && L1<=0 && L2<=0) ||
            (L0>=0 && L1>=0 && L2>=0));
}

typedef tPoint3<float> Point3f;
typedef tPoint3<double> Point3d;
#if defined(REAL_IS_FLOAT)
typedef Point3f Point3;
#else
typedef Point3d Point3;
#endif

GTB_END_NAMESPACE

#ifndef OUTLINE
#include <gtb/graphics/point3.ipp>
#endif

#endif // GTB_POINT3_INCLUDED
