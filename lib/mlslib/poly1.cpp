
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
#include "poly1.h"

using namespace gtb;

MLSLIB_BEGIN_NAMESPACE

template <class REAL>
Poly1<REAL>::Poly1(const Poly1& rhs) : Poly<REAL>(rhs)
{
    std::copy(rhs.coefficients, rhs.coefficients+this->Poly<REAL>::num_coefficientes(), coefficients);
}

template <class REAL>
void Poly1<REAL>::SetCoefficients(areal* p_coeff)
{
    std::copy(p_coeff, p_coeff+this->Poly<REAL>::num_coefficientes(), coefficients);
}

template <class REAL>
typename Poly1<REAL>::areal Poly1<REAL>::get_coefficient(int idx) const
{
    assert(idx < this->Poly<REAL>::num_coefficientes());
    return coefficients[idx];
}

template <class REAL>
void Poly1<REAL>::Normal0(typename Poly1<REAL>::Vector3& N) const
{
    N[0] = -coefficients[1];
    N[1] = -coefficients[2];
    N[2] = 1;
    N.normalize();
}

template <class REAL>
typename Poly1<REAL>::areal Poly1<REAL>::eval(areal x, areal y) const
{
    areal v = 0;
    int cidx = 0;
    for (int p = 0; p <= degree(); ++p)
    {
        for (int ypwr = 0; ypwr <= p; ++ypwr)
        {
            int xpwr = p - ypwr;
            v += coefficients[cidx] * ipow(x, xpwr) * ipow(y, ypwr);
            ++cidx;
        }
    }
    return v;
}

template <class REAL>
void Poly1<REAL>::Normal(areal x, areal y, typename Poly1<REAL>::Vector3& N) const
{
    typename Poly1<REAL>::Vector3 du(1, 0, 0);
    typename Poly1<REAL>::Vector3 dv(0, 1, 0);

    // areal a = coefficients[0]; Unused
    areal b = coefficients[1];
    areal c = coefficients[2];

    du[2] = b;
    dv[2] = c;

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
template <class REAL>
void Poly1<REAL>::curvatures(areal x, areal y, areal& k1, areal& k2)
{
    printf("Poly1::curvatures undefined\n");
}

/*
 * compute the minimal / maximal point of a polynomials
 */
template <class REAL>
void Poly1<REAL>::minpoint(areal& x, areal& y)
{
    printf("Poly1::minpoint undefined\n");
}

/*
 * Compute the tangent plane to a point
 * return two (normalized) vectors that defines it
 */
template <class REAL>
void Poly1<REAL>::tangent_plane(areal x, areal y, 
				typename Poly1<REAL>::Vector3& tx, typename Poly1<REAL>::Vector3& ty)
{
    // areal a = coefficients[0]; Unused
    areal b = coefficients[1];
    areal c = coefficients[2];

    tx[0] = 1;
    tx[1] = 0;
    tx[2] = b;
    ty[0] = 0;
    ty[1] = 1;
    ty[2] = c;

    tx.normalize();
    ty.normalize();
}

template <class REAL>
int Poly1<REAL>::degree() const
{
    return 1;
}

template <class REAL>
Poly<REAL>* Poly1<REAL>::clone() const
{
    return new Poly1(*this);
}

template <class REAL>
void Poly1<REAL>::assign(const Poly<REAL>* rhs)
{
    const Poly1* p1 = static_cast<const Poly1*>(rhs);
    std::copy(p1->coefficients, p1->coefficients+this->Poly<REAL>::num_coefficientes(), coefficients);
}

template<>
void Poly1<double>::write(FILE* f) const
{
    // int NC = num_coefficientes(); Unused
    for (int i = 0; i < 3; ++i) // Should this be i<NC?
    {
        write_double(coefficients[i], f);
    }
}

template<>
void Poly1<float>::write(FILE* f) const
{
    // int NC = num_coefficientes(); Unused
    for (int i = 0; i < 3; ++i) // Should this be i<NC?
    {
        write_float(coefficients[i], f);
    }
}

template<>
void Poly1<double>::read(FILE* f)
{
    // int NC = num_coefficientes(); Unused
    for (int i = 0; i < 3; ++i) // Should this be i<NC?
    {
        double vi;
        read_double(&vi, f);
        coefficients[i] = vi;
    }
}

template<>
void Poly1<float>::read(FILE* f)
{
    // int NC = num_coefficientes(); Unused
    for (int i = 0; i < 3; ++i) // Should this be i<NC?
    {
        float vi;
        read_float(&vi, f);
        coefficients[i] = vi;
    }
}

template class Poly1<float>;
template class Poly1<double>;

MLSLIB_END_NAMESPACE
