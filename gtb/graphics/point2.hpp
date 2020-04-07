
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


#ifndef GTB_POINT2_INCLUDED
#define GTB_POINT2_INCLUDED

#include <gtb/common.hpp>
#include <gtb/real/real.hpp>
#include <gtb/math/math.hpp>
#include <gtb/graphics/vector2.hpp>
#include <vector>

GTB_BEGIN_NAMESPACE

template <class T> class tVector2;

template <class T>
class tPoint2 {
public:
    typedef T value_type;

    tPoint2();
    template <class T2> explicit tPoint2(const tPoint2<T2> &p);
    tPoint2(T px, T py);
    tPoint2(const T p[2]);
    tPoint2 &operator=(const tPoint2 &p);

    T x() const;
    T y() const;

    void set_x(T px);
    void set_y(T py);
    tPoint2 &reset(T px, T py);

    T operator[](unsigned i) const;
    T &operator[](unsigned i);

    bool operator==(const tPoint2 &p) const;
    bool operator!=(const tPoint2 &p) const;
    bool is_zero() const;

    tPoint2 &operator+=(const tVector2<T> &v);
    tPoint2 &operator-=(const tVector2<T> &v);

#if 0
    friend tPoint2 operator+(const tPoint2 &p, const tVector2<T> &v);
    friend tVector2<T> operator-(const tPoint2 &p, const tPoint2 &q);
    friend tPoint2 operator-(const tPoint2 &p, const tVector2<T> &v);
#endif

    tPoint2 &scale(const tPoint2 &origin, T s);
    tPoint2 &translate(const tVector2<T> &v);
    tPoint2 &translate(T dx, T dy);
    tPoint2 &rotate(T theta);
    tPoint2 &rotate(const tPoint2 &origin, T theta);

    static T distance(const tPoint2 &p, const tPoint2 &q);
    static T squared_distance(const tPoint2 &p, const tPoint2 &q);

    static tPoint2 midpoint(const tPoint2 &A, const tPoint2 &B);
    static tPoint2 centroid(const std::vector<tPoint2> &v);
    static tPoint2 inbetween(const tPoint2 &A, const tPoint2 &B,
			     value_type alpha);

    static bool collinear(const tPoint2 &A,
			  const tPoint2 &B,
			  const tPoint2 &C);
    static tVector2<T> normal(const tPoint2 &A,
			      const tPoint2 &B);

    void load() const;
//    void load();

    void read(FILE *fp);
    void write(FILE *fp) const;

    void add(const tPoint2& rhs);
    void multiply(T value);

    friend class tVector2<T>;

#if 0
    friend std::istream &operator>>(std::istream &s, tPoint2 &p);
    friend std::ostream &operator<<(std::ostream &s, const tPoint2 &p);
#endif

    static const tPoint2 ZERO;

protected:
    T _p[2];
};

typedef tPoint2<float> Point2f;
typedef tPoint2<double> Point2d;
#if defined(REAL_IS_FLOAT)
typedef Point2f Point2;// bward compat.
#else
typedef Point2d Point2;// bward compat.
#endif


/*
 * Return > 0 if left-turn
 *        < 0 if right-turn
 *        = 0 if colinear
 * value = area of parallelogram that is defined by the 3 points
 */
template <class T>
T LeftTurn(const tPoint2<T>& p0, const tPoint2<T>& p1, const tPoint2<T>& p2);


/// by Olga - this function returns true if the given
/// point is inside a triangle defined by p0,p1,p2
/// Pre-condition: both point and *this are on xy plane
template <class T>
bool ContainsPoint(const tPoint2<T>& p0, const tPoint2<T>& p1, const tPoint2<T>& p2, const tPoint2<T>& point);


/*
 * Compute the center of the circle defined by the three
 * points p1,p2,p3
 */
template <class T>
void circumcircle(const tPoint2<T>& p1, const tPoint2<T>& p2, const tPoint2<T>& p3, tPoint2<T>& center);

GTB_END_NAMESPACE

#ifndef OUTLINE
#include <gtb/graphics/point2.ipp>
#endif

#endif // GTB_POINT2_INCLUDED
