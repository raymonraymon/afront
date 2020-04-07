
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


#ifndef GTB_VECTOR3_INCLUDED
#define GTB_VECTOR3_INCLUDED

#include <gtb/graphics/point3.hpp>

#include <gtb/common.hpp>
#include <gtb/real/real.hpp>
#include <gtb/math/math.hpp>

GTB_BEGIN_NAMESPACE

template<class T>
class tVector3 {
public:
    typedef T value_type;

    tVector3();
    template <class T2> explicit tVector3(const tVector3<T2> &v);
    explicit tVector3(const tPoint3<T> &p);
    tVector3(value_type vx, value_type vy, value_type vz);
    explicit tVector3(const value_type v[3]);
    explicit tVector3(const value_type v);
    tVector3 &operator=(const tVector3 &v);

    value_type x() const;
    value_type y() const;
    value_type z() const;

    value_type get_x() const;
    value_type get_y() const;
    value_type get_z() const;

    void set_x(value_type vx);
    void set_y(value_type vy);
    void set_z(value_type vz);
    tVector3 &reset(value_type vx, value_type vy, value_type vz);

    value_type operator[](unsigned i) const;
    value_type &operator[](unsigned i);

    bool operator==(const tVector3 &v) const;
    bool operator!=(const tVector3 &v) const;
    bool is_zero() const;

    tVector3 operator-() const;
    tVector3 &operator*=(const tMatrix4<T> &m);
    tVector3 &operator+=(const tVector3 &v);
    tVector3 &operator-=(const tVector3 &v);
    tVector3 &operator*=(value_type a);
    tVector3 &operator/=(value_type a);

    value_type operator*(const tVector3 &v) const;
    value_type operator*(const tPoint3<T> &p) const;

#if 0
    friend tVector3 operator+(const tVector3 &u, const tVector3 &v);
    friend tPoint3<T> operator+(const tPoint3<T> &p, const tVector3 &v);
    friend tVector3 operator-(const tVector3 &u, const tVector3 &v);
    friend tPoint3<T> operator-(const tPoint3<T> &p, const tVector3 &v);
    friend tVector3 operator*(value_type a, const tVector3 &v);
    friend tVector3 operator*(const tVector3 &v, value_type a);
    friend tVector3 operator*(const tMatrix4<T> &m, const tVector3 &v);
    friend tVector3 operator/(const tVector3 &v, value_type a);
#endif

    void element_wise_multiply(const tVector3& rhs);

    value_type dot(const tVector3 &v) const;
    value_type dot(const tPoint3<T> &p) const;
    tVector3 cross(const tVector3 &v) const;
    value_type length() const;
    value_type squared_length() const;
    tVector3 &normalize();
    tVector3 normalized() const;
    bool is_normalized() const;
    tVector3 &scale(value_type a);
    tVector3 &rotate(const tVector3 &axis, value_type theta);
    tVector3 &transform(const tMatrix4<T> &m);
    tVector3 &affine_transform(const tMatrix4<T> &m);
    void flip();
    tVector3 flipped() const;

    void load_as_normal() const;

    void read(FILE *fp);
    void write(FILE *fp) const;

    static tMatrix4<T> rotation(const tVector3 &from, const tVector3 &to);

#if 0
    friend class tPoint3<T>;
    friend void swap(tVector3 &a, tVector3 &b);
    friend std::istream &operator>>(std::istream &s, tVector3 &v);
    friend std::ostream &operator<<(std::ostream &s, const tVector3 &v);
#endif

    static const tVector3 VECTOR3_ZERO;
    static const tVector3 VECTOR3_POSITIVE_X;
    static const tVector3 VECTOR3_NEGATIVE_X;
    static const tVector3 VECTOR3_POSITIVE_Y;
    static const tVector3 VECTOR3_NEGATIVE_Y;
    static const tVector3 VECTOR3_POSITIVE_Z;
    static const tVector3 VECTOR3_NEGATIVE_Z;

protected:
    value_type _v[3];
};

typedef tVector3<float> Vector3f;
typedef tVector3<double> Vector3d;
#if defined(REAL_IS_FLOAT)
typedef Vector3f Vector3;
#else
typedef Vector3d Vector3;
#endif

GTB_END_NAMESPACE

#ifndef OUTLINE
#include <gtb/graphics/vector3.ipp>
#endif

#endif // GTB_VECTOR3_INCLUDED
