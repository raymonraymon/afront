
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


#ifndef __POLY_H
#define __POLY_H

#include "mlslibdefs.h"
#include <gtb/gtb.hpp>

MLSLIB_BEGIN_NAMESPACE

/*
 * generic polynomial class
 */
template <class REAL>
class Poly
{
public:
    typedef REAL areal;
    typedef gtb::tVector3<REAL> Vector3;

    virtual ~Poly() {}

    virtual void SetCoefficients(areal* p_coeff)=0;
    virtual areal eval0() const=0;
    virtual void Normal0(Vector3& N) const=0;

    virtual Poly* clone() const = 0;
    virtual void assign(const Poly* rhs) = 0;
    virtual areal get_coefficient(int idx) const = 0;

    virtual areal eval(areal x, areal y) const = 0;
    virtual void Normal(areal x, areal y, Vector3& N) const=0;
    virtual void curvatures(areal x, areal y, areal& k1, areal& k2) {};
    virtual void minpoint(areal& x, areal& y){};
    virtual void tangent_plane(areal x, areal y, Vector3& tx, Vector3& ty){};

    virtual int degree() const=0;
    virtual int num_coefficientes() const { return (degree() + 1)*(degree() + 2)/2; }

    virtual void write(FILE* f) const = 0;
    virtual void read(FILE* f) = 0;
};

template <class REAL>
void write_poly(FILE* f, const Poly<REAL>* c);

template <class REAL>
Poly<REAL>* read_poly(FILE* f);

MLSLIB_END_NAMESPACE

#endif // __POLY_H
