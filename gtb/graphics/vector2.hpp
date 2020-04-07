
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


#ifndef GTB_VECTOR2_INCLUDED
#define GTB_VECTOR2_INCLUDED

#include <gtb/common.hpp>
#include <gtb/real/real.hpp>
#include <gtb/math/math.hpp>
#include <gtb/graphics/point2.hpp>

GTB_BEGIN_NAMESPACE


template <class T> class tPoint2;

//class tPoint2<T>;
//class Vector3;

template <class T>
class tVector2 {
public:
    typedef T value_type;

    tVector2();
    tVector2(const tVector2 &v);
    tVector2(T vx, T vy);
    tVector2(const T v[2]);
    tVector2 &operator=(const tVector2 &v);
    template <class T2>
    tVector2 (const tVector2<T2> &v) {
	_v[0] = T(v[0]);
	_v[1] = T(v[1]);
    }

    T x() const;
    T y() const;
    T get_x() const;
    T get_y() const;

    void set_x(T vx);
    void set_y(T vy);
    tVector2 &reset(T vx, T vy);

    T operator[](unsigned i) const;
    T &operator[](unsigned i);

    bool operator==(const tVector2 &v) const;
    bool operator!=(const tVector2 &v) const;
    bool is_zero() const;

    tVector2 operator-() const;
    tVector2 &operator+=(const tVector2 &v);
    tVector2 &operator-=(const tVector2 &v);
    tVector2 &operator*=(T a);
    tVector2 &operator/=(T a);

    value_type operator*(const tVector2 &v) const;
    value_type operator*(const tPoint2<T> &v) const;

#if 0
    friend tVector2 operator+(const tVector2 &u, const tVector2 &v);
    friend tPoint2<T> operator+(const tPoint2<T> &p, const tVector2 &v);
    friend tVector2 operator-(const tVector2 &u, const tVector2 &v);
    friend tPoint2<T> operator-(const tPoint2<T> &p, const tVector2 &v);
    friend tVector2 operator*(T a, const tVector2 &v);
    friend tVector2 operator*(const tVector2 &v, T a);
    friend tVector2 operator/(const tVector2 &v, T a);
#endif

    T dot(const tVector2 &v) const;
    T dot(const tPoint2<T> &v) const;
    tVector3<T> cross(const tVector2 &v) const;
    T length() const;
    T squared_length() const;
    tVector2 &normalize();
    bool is_normalized() const;
    tVector2 &scale(T a);
    tVector2 &rotate(T theta);
    tVector2 normal() const;

    friend class tPoint2<T>;
#if 0
    friend std::istream &operator>>(std::istream &s, tVector2 &v);
    friend std::ostream &operator<<(std::ostream &s, const tVector2 &v);
#endif

    static tVector2 ZERO;
    static tVector2 POSITIVE_X;
    static tVector2 NEGATIVE_X;
    static tVector2 POSITIVE_Y;
    static tVector2 NEGATIVE_Y;

    void read(FILE *fp) {
	type_traits<T>::read(&_v[0], fp);
	type_traits<T>::read(&_v[1], fp);
    }
    void write(FILE *fp) const {
	type_traits<T>::write(_v[0], fp);
	type_traits<T>::write(_v[1], fp);
    }

protected:
    T _v[2];
};

GTB_GENERATE_CLASS_TYPEDEFS(Vector2)

GTB_END_NAMESPACE

#ifndef OUTLINE
#include <gtb/graphics/vector2.ipp>
#endif

#endif // GTB_VECTOR2_INCLUDED
