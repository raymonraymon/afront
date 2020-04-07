
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


#ifndef __TMAT3_HPP
#define __TMAT3_HPP

#include <iostream>
#include <gtb/graphics/vector3.hpp>
#include <gtb/error/error.hpp>

GTB_BEGIN_NAMESPACE

template <class T>
class tmat3
{
protected:

    tVector3<T> v[3];

public:
    typedef T value_type;
    typedef tVector3<T> Vector3;

    // Constructors

    tmat3();
    tmat3(const tVector3<T>& v0, const tVector3<T>& v1, const tVector3<T>& v2);
    tmat3(const double d);
    tmat3(const tmat3& m);

    // Assignment operators

    tmat3& operator	= ( const tmat3& m );	    // assignment of a tmat3
    tmat3& operator += ( const tmat3& m );	    // incrementation by a tmat3
    tmat3& operator -= ( const tmat3& m );	    // decrementation by a tmat3
    tmat3& operator *= ( const double d );	    // multiplication by a constant
    tmat3& operator /= ( const double d );	    // division by a constant
    tVector3<T>& operator [] ( int i);		    // indexing
    const tVector3<T>& operator [] ( int i) const;		    // indexing

    // special functions

    tmat3 transpose() const;			    // transpose
    tmat3 inverse();				    // inverse
    //tmat3& apply(V_FCT_PTR fct);		    // apply a func. to each element

#if 0
    // friends
    friend tmat3 operator - (const tmat3& a);			    // -m1
    friend tmat3 operator + (const tmat3& a, const tmat3& b);	    // m1 + m2
    friend tmat3 operator - (const tmat3& a, const tmat3& b);	    // m1 - m2
    friend tmat3 operator * (const tmat3& a, const tmat3& b);	    // m1 * m2
    friend tmat3 operator * (const tmat3& a, const double d);	    // m1 * 3.0
    friend tmat3 operator * (const double d, const tmat3& a);	    // 3.0 * m1
    friend tmat3 operator / (const tmat3& a, const double d);	    // m1 / 3.0
    friend int operator == (const tmat3& a, const tmat3& b);	    // m1 == m2 ?
    friend int operator != (const tmat3& a, const tmat3& b);	    // m1 != m2 ?
    //#ifndef WINDOWS
    friend std::ostream& operator << (std::ostream& s, tmat3& m);	    // output to stream
    friend std::istream& operator >> (std::istream& s, tmat3& m);	    // input from strm.
    //#endif
    friend void swap(tmat3& a, tmat3& b);			    // swap m1 & m2
#endif

    value_type determinant() const;
    value_type MINOR(int r0, int r1, int c0, int c1) const;
    bool IsSingular() const;

    //
    // From WildMagic
    //
    void FromEulerAnglesXYZ(value_type fYAngle, value_type fPAngle, value_type fRAngle);
    bool ToEulerAnglesXYZ(value_type& rfXAngle, value_type& rfYAngle, value_type& rfZAngle) const;


    // necessary friend declarations

#if 0
    friend tVector3<T> operator * (const tmat3& a, const tVector3<T>& v);	    // linear transform
    friend Point3 operator * (const tmat3& a, const Point3& v);	    // linear transform
#endif
};

typedef tmat3<float> mat3f;
typedef tmat3<double> mat3d;
#if defined(REAL_IS_FLOAT)
typedef mat3f mat3;
#else
typedef mat3d mat3;
#endif

//extern tmat3 tmat3I;
//extern tmat3 tmat3Zero;

#include "mat3.inl"

GTB_END_NAMESPACE

#endif // __TMAT3_HPP
