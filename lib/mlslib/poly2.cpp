
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


#include "stdafx.h"
#include "mlslibdefs.h"
#include "poly2.h"

using namespace gtb;

MLSLIB_BEGIN_NAMESPACE

template<class REAL>
Poly2<REAL>::Poly2(const Poly2& rhs) : Poly<REAL>(rhs)
{
    std::copy(rhs.coefficients, rhs.coefficients+this->Poly<REAL>::num_coefficientes(), coefficients);
}

template<class REAL>
void Poly2<REAL>::SetCoefficients(areal* p_coeff)
{
    std::copy(p_coeff, p_coeff+this->Poly<REAL>::num_coefficientes(), coefficients);
}

template<class REAL>
typename Poly2<REAL>::areal Poly2<REAL>::get_coefficient(int idx) const
{
    assert(idx < this->Poly<REAL>::num_coefficientes());
    return coefficients[idx];
}

template<class REAL>
void Poly2<REAL>::Normal0(Vector3& N) const
{
    N[0] = -coefficients[1];
    N[1] = -coefficients[2];
    N[2] = 1;
    N.normalize();
}

template<class REAL>
typename Poly2<REAL>::areal Poly2<REAL>::eval(areal x, areal y) const
{
    areal v = 0;
    int cidx = 0;
//    int deg = degree();
    static const int deg=2;
    for (int p = 0; p <= deg; ++p)
    {
        areal ypy=1.0;
        for (int ypwr = 0; ypwr <= p; ++ypwr)
        {
            int xpwr = p - ypwr;
            v += coefficients[cidx] * ipow(x, xpwr) * ypy;
            ypy *= y;
            ++cidx;
        }
    }
    return v;
}

template<class REAL>
void Poly2<REAL>::Normal(areal x, areal y, Vector3& N) const
{
    Vector3 du(1, 0, 0);
    Vector3 dv(0, 1, 0);

    //    areal a = coefficients[0]; Unused
    areal b = coefficients[1];
    areal c = coefficients[2];
    areal d = coefficients[3];
    areal e = coefficients[4];
    areal f = coefficients[5];

    du[2] = b + 2*d*x + e*y;
    dv[2] = c + e*x + 2*f*y;

    N = du.cross(dv);
    N.normalize();
}

/*
 * Compute the minimal, maximal curvature at x,y
 *
 * return:
 *  k1 - minimal curvature
 *  k2 - maximal curvature
 */
template<class REAL>
void Poly2<REAL>::curvatures(areal x, areal y, areal& k1, areal& k2)
{
    //    areal a = coefficients[0]; Unused
    areal b = coefficients[1];
    areal c = coefficients[2];
    areal d = coefficients[3];
    areal e = coefficients[4];
    areal f = coefficients[5];

    areal nx = b + 2*d*x + e*y;
    areal ny = c + e*x + 2*f*y;
    areal nlen = sqrt(nx*nx+ny*ny+1);

    areal b11 = 2*d/nlen;
    areal b12 = e/nlen;
    areal b22 = 2*f/nlen;
    areal disc = (b11+b22)*(b11+b22) - 4*(b11*b22-b12*b12);
    assert (disc>=0);
    areal disc2 = sqrt(disc);
    k1 = (b11+b22+disc2)/2.0;
    k2 = (b11+b22-disc2)/2.0;
    if (gtb::absv(k1) > gtb::absv(k2)) std::swap(k1, k2);
}

/*
 * compute the minimal / maximal point of a polynomials
 */
template<class REAL>
void Poly2<REAL>::minpoint(areal& x, areal& y)
{
    //    areal a = coefficients[0]; Unused
    areal b = coefficients[1];
    areal c = coefficients[2];
    areal d = coefficients[3];
    areal e = coefficients[4];
    areal f = coefficients[5];

    areal d2 = d*d;
    areal e2 = e*e;

    areal quatx = 8*f*d2-2*d*e2;
    if (gtb::absv(quatx) < 1e-8) x = 0;
    else x = (b*e2-4*b*d*f+b*e2-2*d*c*e)/quatx;
    areal quaty = 4*f*d-e2;
    if (gtb::absv(quaty) < 1e-8) y = 0;
    else y = (b*e-2*d*c)/quaty;
}

/*
 * Compute the tangent plane to a point
 * return two (normalized) vectors that defines it
 */
template<class REAL>
void Poly2<REAL>::tangent_plane(areal x, areal y, Vector3& tx, Vector3& ty)
{
    //    areal a = coefficients[0]; Unused
    areal b = coefficients[1];
    areal c = coefficients[2];
    areal d = coefficients[3];
    areal e = coefficients[4];
    areal f = coefficients[5];

    tx[0] = 1;
    tx[1] = 0;
    tx[2] = b + 2*d*x + e*y;
    ty[0] = 0;
    ty[1] = 1;
    ty[2] = c + e*x + 2*f*y;

    tx.normalize();
    ty.normalize();
}

template<class REAL>
int Poly2<REAL>::degree() const
{
    return 2;
}

template<class REAL>
Poly<REAL>* Poly2<REAL>::clone() const
{
    return new Poly2(*this);
}

template<class REAL>
void Poly2<REAL>::assign(const Poly<REAL>* rhs)
{
    const Poly2* p2 = static_cast<const Poly2*>(rhs);
    std::copy(p2->coefficients, p2->coefficients+this->Poly<REAL>::num_coefficientes(), coefficients);
}

template<>
void Poly2<double>::write(FILE* f) const
{
    // int NC = num_coefficientes(); Unused
    for (int i = 0; i < 6; ++i) // should this be i<NC?
    {
        write_double(coefficients[i], f);
    }
}

template<>
void Poly2<float>::write(FILE* f) const
{
    // int NC = num_coefficientes(); Unused
    for (int i = 0; i < 6; ++i) // should this be i<NC?
    {
        write_float(coefficients[i], f);
    }
}

template<>
void Poly2<double>::read(FILE* f)
{
    // int NC = num_coefficientes(); Unused
    for (int i = 0; i < 6; ++i) // should this be i<NC?
    {
        double vi;
        read_double(&vi, f);
        coefficients[i] = vi;
    }
}

template<>
void Poly2<float>::read(FILE* f)
{
    // int NC = num_coefficientes(); Unused
    for (int i = 0; i < 6; ++i) // should this be i<NC?
    {
        float vi;
        read_float(&vi, f);
        coefficients[i] = vi;
    }
}

/*
 * Compute the curvature at the minimal point of the poly
 */
template<class REAL>
void Poly2<REAL>::curvatures_at_minimum(areal& k1, areal& k2)
{
    // find the minimal point
    areal b = coefficients[1];
    areal c = coefficients[2];
    areal d = coefficients[3];
    areal e = coefficients[4];
    areal f = coefficients[5];
    areal det = 4*d*f-e*e;
    if (absv(det) < 1e-10)
    {
        k1 = k2 = 0;
    }
    else
    {
        areal x = (-2*f*b + c*e) / det;
        areal y = (b*e - 2*d*e) / det;
        curvatures(x,y,k1,k2);
    }
}

template class Poly2<float>;
template class Poly2<double>;

MLSLIB_END_NAMESPACE
