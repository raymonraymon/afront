
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


#ifndef __GTB_POLY_H
#define __GTB_POLY_H

#include <gtb/real/real.hpp>
#include <gtb/math/amat.h>
#include <gtb/math/mathtools.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include <gtb/mp/mp.h>
#include <gtb/mp/mp_tuple.h>
#include <iomanip>
#include <gtb/math/linearsolvers.h>

GTB_BEGIN_NAMESPACE

enum poly_lls_solver {
	POLY_LLS_SVD,
	POLY_LLS_CHOLESKY,
	POLY_LLS_HOUSEHOLDER,
};


// This is just a 1D polynomial on an arbitrary basis and degree
template <int Degree,
          class Real,
	  template <typename> class Basis>
class Poly1: public Basis<Real> 
{
public:
    Poly1() {}
    Real Eval(Real x) const;
    void Set(int ix, Real x) {
	_coefs[imap(ix)] = x;
    }
    Real Get(int ix) {
	return _coefs[imap(ix)];
    }
    void Clear() {
	std::fill(_coefs, _coefs + CoefCount, 0);
    }

    void LeastSquares(const std::vector<Real> &xs,
		      const std::vector<Real> &fxs);

    void write(std::ostream &os);

protected:
    enum { CoefCount = (Degree+1) };
    Real _coefs[CoefCount];
    int imap(int ix) const { return ix; }
};

// This is just a 2D polynomial on an arbitrary basis and degree
template <int Degree,
          class Real,
	  template <typename> class Basis>
class Poly2: public Basis<Real> 
{
public:
    Poly2() {}
    Real Eval(Real x, Real y) const;
    void Set(int ix, int iy, Real x) {
	_coefs[imap(ix, iy)] = x;
    }
    Real Get(int ix, int iy) {
	return _coefs[imap(ix, iy)];
    }
    void Clear() {
	std::fill(_coefs, _coefs + CoefCount, 0);
    }

	void LeastSquares(const std::vector<Real> &xs, const std::vector<Real> &ys,
		      const std::vector<Real> &fxs);

protected:
    enum { CoefCount = (Degree+1)*(Degree+2)/2 };
    Real _coefs[CoefCount];
    int imap(int ix, int iy) const { 
	iy = ix+iy;
	return ix + (iy+1)*iy/2; 
    }
};

// This is just a 3D polynomial on an arbitrary basis and degree
template <int Degree,
	  class Real,
	  template <typename> class Basis>
class Poly3: public Basis<Real>
{
public:
    Poly3() {}
    Real Eval(Real x, Real y, Real z) const;
    void Set(int ix, int iy, int iz, Real x) {
	_coefs[imap(ix, iy, iz)] = x;
    }
    Real Get(int ix, int iy, int iz) {
	return _coefs[imap(ix, iy, iz)];
    }
    void Clear() {
	std::fill(_coefs, _coefs + CoefCount, 0);
    }

	void LeastSquares(const std::vector<Real> &xs, const std::vector<Real> &ys, const std::vector<Real> &zs,
		      const std::vector<Real> &fxs);

	void LeastSquares(const std::vector<Real> &xs, const std::vector<Real> &ys, const std::vector<Real> &zs,
		      const std::vector<Real> &fxs, const std::vector<Real> &weights, bool normeqns, poly_lls_solver solver);

    //protected:
    enum { CoefCount = (Degree+1)*(Degree+2)*(Degree+3)/6 };
    Real _coefs[CoefCount];
    int imap(int ix, int iy, int iz) const {
	iy = ix+iy;
	iz = iy+iz;
	return ix + (iy+1)*iy/2 + (iz+2)*(iz+1)*iz/6; 
    }
};

template <class Real>
class Monomial {
public:
    // Evaluate single term
    static Real TermEval(int n, Real v)
    {
	return ipow(v, n);
    };
};

template <class Real>
class Legendre {
public:
    static Real TermEval(int n, Real x) {
	Real result = (Real)0;
	for (int k=0; k<=n/2; ++k)
	    result += ipow(-1, k) 
		* fchoose(n, k)
		* fchoose(2*n-2*k, n)
		* ipow(x, n-2*k);
	return ipow((Real)0.5, n) * result;
    };
};

template <int Dimension, int Degree, class Real = real_t,
	  template <typename> class Basis = Monomial>
class Poly;

#include "poly_trimapper.inl"

#include "poly_termeval.inl"

template 
    <int Dimension, 
     int Degree, 
     class Real,
     template <typename> class Basis>
    class Poly
{
public:
    enum { CoefCount = mp::Choose<Degree+Dimension, Degree>::Value };
    void Clear() {
	fill(_coefs, _coefs + CoefCount, (Real)0);
    }

    Real Eval(const typename mp::RealTuple<Dimension>::Type &x) const 
    {
	return 
	    PolyEval<Dimension, Degree, Real, Basis, Degree, 0>::Eval
	    (*this, x, gtb::mp::Tuple<int, 0>());
    }

    Real Eval(Real x) const
    {
	return Eval(GTB_MP_MAKE_REAL_TUPLE_1(x));
    }
    Real Eval(Real x, Real y) const
    {
	return Eval(GTB_MP_MAKE_REAL_TUPLE_2(x,y));
    }
    Real Eval(Real x, Real y, Real z) const
    {
	return Eval(GTB_MP_MAKE_REAL_TUPLE_3(x,y,z));
    }

    Real Get(const typename mp::IntTuple<Dimension>::Type &x) const
    {
	return _coefs[imap(x)];
    }

    void Set(const typename mp::IntTuple<Dimension>::Type &x,
	     Real v) 
    {
	_coefs[imap(x)] = v;
    }

    void Set(int a, Real t) 
    { Set(GTB_MP_MAKE_INT_TUPLE_1(a), t); }

    Real Get(int a) const 
    { return Get(GTB_MP_MAKE_INT_TUPLE_1(a)); }

    void Set(int a, int b, Real t) 
    { Set(GTB_MP_MAKE_INT_TUPLE_2(a,b), t); }

    Real Get(int a, int b) const 
    { return Get(GTB_MP_MAKE_INT_TUPLE_2(a,b)); }

    void Set(int a, int b, int c, Real t) 
    { Set(GTB_MP_MAKE_INT_TUPLE_3(a,b,c), t); }

    Real Get(int a, int b, int c) const 
    { return Get(GTB_MP_MAKE_INT_TUPLE_3(a,b,c)); }

    Real _coefs[CoefCount];

    void LeastSquares(const gtb::AMat<Real> &xs,
		      const gtb::AVec<Real> &fxs);

    void LeastSquares(const gtb::AMat<Real> &xs,
		      const gtb::AVec<Real> &fxs,
		      const gtb::AVec<Real> &weights);

    void write(std::ostream &os) const;
    int imap(const mp::Tuple<int, Dimension> &x) const {
	return triangular_map(x);
    };
    int imap(int *x) const {
	return triangular_map<Dimension>(x);
    };
};

#include "poly.inl"

GTB_END_NAMESPACE

#endif // __POLY_H
