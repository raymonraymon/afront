
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


#ifndef __AMAT_H
#define __AMAT_H

/*
 * A general purpose matrix and vector classes
 * incomplete, but satisfies my needs
 */

#include <gtb/common.hpp>
#include <cassert>

GTB_BEGIN_NAMESPACE

template<typename T>
class AVec
{
public:
    typedef T value_type;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef pointer iterator;
    typedef const_pointer const_iterator;

    AVec(int n, const_pointer p=0);
    AVec(const AVec& rhs);
    ~AVec();

    reference operator()(int n);
    const_reference operator()(int n) const;

    reference operator[](int n);
    const_reference operator[](int n) const;

    AVec& operator=(const AVec& rhs);

    void assign(const pointer rhs);

    int size() const { return N; }
    void resize(int n);

    void set(const T& value);

    operator pointer() { return v; }
    operator const_pointer() const { return v;}

    iterator begin() { return v; }
    const_iterator begin() const { return v; } 
    iterator end() { return v+N; }
    const_iterator end() const { return v+N; }

    pointer fortran_pointer() { return v-1; }

    T dot(const AVec& rhs) const;

protected:
    pointer v;
    int N;

    void allocate();
    void cleanup();
    void copy(const_pointer rhs);
};

template <class T>
class AMat
{
public:
    typedef T value_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef value_type* pointer;
    typedef const value_type*  const_pointer;
    typedef pointer* pointer2;
    typedef const  const_pointer* const_pointer2;

    AMat(int m, int n, const_pointer v, bool ftn=false);
    AMat(int m, int n, const_pointer2 A=0);
    AMat(const AMat& rhs);
    ~AMat();

    AMat operator *(const AMat& rhs) const;
    AVec<T> operator *(const AVec<T>& rhs) const;
    void transpose();
    AMat transposed() const;

    reference operator() (int i, int j);
    const_reference operator() (int i, int j) const;

    AMat& operator=(const AMat& rhs);

    void assign(const_pointer v, bool ftn=false);
    void assign(const_pointer v, int m, int n);

    void assign(const_pointer2 A);
    void assign(const_pointer2 A, int m, int n);

    void assign1(const_pointer2 A);
    void assign1(const_pointer2 A, int m, int n);

    void assign_diag(const AVec<T>& v);

    int rows() const { return M; }
    int columns() const { return N; }
    void resize(int m, int n);

    void set(const T& value);

    void get_row(int r, AVec<T>& vrow);
    void get_column(int c, AVec<T>& vcol);

protected:
    value_type** m;
    int M,N;

    void allocate();
    void cleanup();
    void copy(const_pointer2 A);
};    


/*
 * Utility I/O functions
 */
template <class T>
std::ostream& operator << (std::ostream& s, const AVec<T>& v);
template <class T>
std::ostream& operator << (std::ostream& s, const AMat<T>& m);

template <class T>
AVec<T> flatten(const AMat<T>& m);

template <class T>
AVec<T> fortran_flatten(const AMat<T>& m);

/*
 * Transform a matrix to a flat vector
 * transpose
 */
template <class T>
inline AVec<T> fortran_flatten(const AMat<T>& m)
{
    AVec<T> r(m.rows() * m.columns());
    int p=0;
	for (int j = 0; j < m.columns(); ++j)
	{
	     for (int i = 0; i < m.rows(); ++i, ++p)
		 {
			 r(p) = m(i,j);
		 }
	}

    return r;
}

//
// Input: matrix flattned by fortran_flatten (i.e. row-major)
// output: a matrix
//
// Assume M is in the correct size
//
template <class T>
void fortran_expand(const AVec<T>& fM, AMat<T>& M);

template <class T>
void diag(const AMat<T>& m, AVec<T> v);


/*
 * Commonly used
 */
typedef AVec<float> AVecf;
typedef AVec<double> AVecd;
typedef AMat<float> AMatf;
typedef AMat<double> AMatd;

/*
 * Return x'*M*x
 */
template <class T>
T xtMx(const AVec<T>& x, const AMat<T>& M);

/*
 * Compute M'M
 */
template <class T>
AMat<T> MtM(const AMat<T>& M);

/*
 * Linear algebra tools
 */
void SVD(
    const AMatd& A,     // matrix to solve
    AMatd& U,           // result U
    AVecd& W,           // result diagonal "matrix" W
    AMatd& Vt           // result V^t
    );

AMatd SVDInvert(const AMatd& A);

void SVDSolve(
    const AMatd& A,     // matrix to solve
    const AVecd& b,
    AVecd& x            // Result
    );

//
// Compute the eigen vectors and values
// of a square matrix
//
void EigenUpperTriangular(const AMatd& A, AMatd& V, AVecd& evalues);

#include "amat.inl"
GTB_END_NAMESPACE

#endif //__AMAT_H
