
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


#include <gtb/gtb.hpp>
#include <assert.h>
#include <gtb/math/amat.h>

//#define USEMKL

#ifdef USENR
#include "NR/nr.h"
#else
#ifdef USEMKL
#include <mkl.h>
#elif __APPLE__
//#include <f2c.h>
#include <vecLib/clapack.h>
#else
#include <f2c.h>
extern "C" {
#include <clapack.h>
}
#endif
#endif

#include "amat.h"

GTB_BEGIN_NAMESPACE

/*------------------ AVec -----------------------*/
template<class T>
AVec<T>::AVec(int n, const_pointer p) : N(n)
{
    allocate();
    if (p) copy(p);
}


template<class T>
void AVec<T>::assign(const pointer rhs)
{
    copy(rhs);
}

template<class T>
void AVec<T>::set(const T& value)
{
    for (int i = 0; i < N; ++i)
    {
        v[i] = value;
    }
}

template<class T>
void AVec<T>::allocate()
{
    v = new value_type[N];
}

template<class T>
void AVec<T>::copy(typename AVec<T>::const_pointer rhs)
{
    std::copy(rhs, rhs+N, v);
}

template<class T>
AVec<T>& AVec<T>::operator=(const AVec<T>& rhs)
{
    resize(rhs.size());
    copy(rhs);
    return *this;
}

template<class T>
void AVec<T>::resize(int n)
{
    if (N != n)
    {
        cleanup();
        N=n;
        allocate();
    }
}

template<class T>
T AVec<T>::dot(const AVec& rhs) const
{
    assert(size() == rhs.size());
    int N = size();
    T result = 0;
    for (int i = 0; i < N; ++i)
    {
        result += (*this)(i) * rhs(i);
    }
    return result;
}

/*------------------ AMat ---------------------*/
template<class T>
AMat<T>::AMat(int row, int col, const_pointer v, bool ftn) : M(row), N(col)
{
    allocate();
    if (v) assign(v, ftn);
}

template<class T>
AMat<T>::AMat(int row, int col, const_pointer2 A) : M(row), N(col)
{
    allocate();
    if (A) copy(A);
}

template<class T>
AMat<T>::AMat(const AMat<T>& rhs) : M(rhs.M), N(rhs.N)
{
    allocate();
    copy(rhs.m);
}

template<class T>
AMat<T>::~AMat()
{
    cleanup();
}

template<class T>
AMat<T> AMat<T>::operator *(const AMat<T>& rhs) const
{
    assert(N == rhs.M);

    AMat<T> r(M, rhs.N);

    for (int i = 0; i < r.M; ++i)
    {
        for (int j = 0; j < r.N; ++j)
        {
            r.m[i][j] = 0;
            for (int k = 0; k < N; ++k)
            {
                r.m[i][j] += m[i][k] * rhs.m[k][j];
            }
        }
    }

    return r;
}

template<class T>
AVec<T> AMat<T>::operator *(const AVec<T>& rhs) const
{
    assert(N == rhs.size());

    AVec<T> r(M);

    for (int i = 0; i < M; ++i)
    {
        r(i) = 0;
        for (int j = 0; j < N; ++j)
        {
            r[i] += m[i][j]*rhs[j];
        }
    }

    return r;
}

template<class T>
void AMat<T>::transpose()
{
    if (M==N)
    {
        for (int i = 0; i < M; ++i)
        {
            for (int j = i; j < N; ++j)
            {
                if (i != j) std::swap(m[i][j], m[j][i]);
            }
        }
    }
    else
    {
        AMat temp = transposed();
        *this = temp;
    }
}


template<class T>
AMat<T> AMat<T>::transposed() const
{
    int N1 = M;
    int M1 = N;
    AMat result(M1, N1);
    for (int i = 0; i < M1; ++i)
        for (int j = 0; j < N1; ++j)
            result(i,j) = (*this)(j,i);
    return result;
}

template<class T>
AMat<T>& AMat<T>::operator=(const AMat& rhs)
{
    cleanup();
    M=rhs.M;
    N=rhs.N;
    allocate();
    copy(rhs.m);
    return *this;
}

/*
 * If ftn==true, transpose the matrix
 * fotran matrices are flat and transposed
 */
template<class T>
void AMat<T>::assign(typename AMat<T>::const_pointer v, bool ftn)
{
    int p=0;

	if (!ftn) {
		for (int i = 0; i < M; ++i)
		{
			for (int j = 0; j < N; ++j, ++p)
			{
				m[i][j] = v[p];
			}
		}
	} else {
		for (int j = 0; j < N; ++j)
		{
			for (int i = 0; i < M; ++i, ++p)
			{
				m[i][j] = v[p];
			}
		}
	}
}

template<class T>
void AMat<T>::assign(typename AMat<T>::const_pointer v, int row, int col)
{
    cleanup();
    M=row;
    N=col;
    allocate();

    assign(v);
}

template<class T>
void AMat<T>::assign(typename AMat<T>::const_pointer2 A)
{
    copy(A);
}

template<class T>
void AMat<T>::assign(typename AMat<T>::const_pointer2 A, int row, int col)
{
    cleanup();
    M=row;
    N=col;
    allocate();

    assign(A);
}

/*
 * Copy from a matrix that is defined in the range of 1..m x 1..n
 */
template<class T>
void AMat<T>::assign1(typename AMat<T>::const_pointer2 A)
{
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < N; ++j)
            m[i][j] = A[i+1][j+1];
}

/*
 * Copy from a matrix that is defined in the range of 1..m x 1..n
 */
template<class T>
void AMat<T>::assign1(typename AMat<T>::const_pointer2 A, int row, int col)
{
    cleanup();
    M=row;
    N=col;
    allocate();

    assign1(A);
}

/*
 * Create a diagonal matrix from a vector
 */
template<class T>
void AMat<T>::assign_diag(const AVec<T>& v)
{
    for (int i = 0; i < M; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            if (i == j) m[i][j] = v[i];
            else m[i][j] = 0;
        }
    }
}

template<class T>
void AMat<T>::allocate()
{
    m = new value_type*[M];
    for (int i = 0; i < M; ++i)
    {
        m[i] = new value_type[N];
    }
}

template<class T>
void AMat<T>::cleanup()
{
    for (int i = 0; i < M; ++i)
    {
        delete[] m[i];
    }
    delete[] m;
}

template<class T>
void AMat<T>::copy(typename AMat<T>::const_pointer2 A)
{

    for (int i = 0; i < M; ++i)
        for (int j = 0; j < N; ++j)
            m[i][j] = A[i][j];
}

template<class T>
void AMat<T>::resize(int row, int col)
{
    if ((row != M) || (col != N))
    {
        cleanup();
        M=row;
        N=col;
        allocate();
    }
}

template<class T>
void AMat<T>::set(const T& value)
{
    for (int i = 0; i < M; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
                        m[i][j] = value;
                }
        }
}


/*
 * Return the rth row as a vector
 */
template <class T>
void AMat<T>::get_row(int r, AVec<T>& vrow)
{
    vrow.resize(columns());
    for (int i = 0; i < columns(); ++i)
    {
        vrow(i) = (*this)(r, i);
    }
}

/*
 * Return the cth column as a vector
 */
template <class T>
void AMat<T>::get_column(int c, AVec<T>& vcol)
{
    vcol.resize(rows());
    for (int i = 0; i < rows(); ++i)
    {
        vcol(i) = (*this)(i, c);
    }
}

/*------------------ Utilities -----------------------*/

template <class T>
std::ostream& operator << (std::ostream& s, const AVec<T>& v)
{
    for (int i = 0; i < v.size(); ++i)
    {
        s << v(i) << ' ';
    }

    return s;
}


template <class T>
std::ostream& operator << (std::ostream& s, const AMat<T>& m)
{
    for (int i = 0; i < m.rows(); ++i)
    {
        for (int j = 0; j < m.columns(); ++j)
            s << m(i,j) << ' ';
        s << '\n';
    }

    return s;
}

/*
 * Transform a matrix to a flat vector
 */
template <class T>
AVec<T> flatten(const AMat<T>& m)
{
    AVec<T> r(m.rows() * m.columns());
    int p=0;
    for (int i = 0; i < m.rows(); ++i)
    {
        for (int j = 0; j < m.columns(); ++j, ++p)
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
void fortran_expand(const AVec<T>& fM, AMat<T>& M)
{
    assert(M.rows() * M.columns() == fM.size());

    int p = 0;
    for (int j = 0; j < M.columns(); ++j)
    {
        for (int i = 0; i < M.rows(); ++i, ++p)
        {
            M(i,j) = fM(p);
        }
    }
}


/*
 * Copy the diagonal of a matrix to a vector, assume the size are correct
 */
template <class T>
void diag(const AMat<T>& m, AVec<T> v)
{
    assert(m.rows() == m.columns() == v.size());
    int N = v.size();
    for (int i = 0; i < N; ++i)
    {
        v(i) = m(i,i);
    }
}

// Instantiate for float and double
template class AVec<float>;
template class AVec<double>;
template class AMat<float>;
template class AMat<double>;

using namespace std;
template ostream& operator << <float> (ostream& s, const AVec<float>& v);
template ostream& operator << <double> (ostream& s, const AVec<double>& v);
template ostream& operator << <float> (ostream& s, const AMat<float>& m);
template ostream& operator << <double> (ostream& s, const AMat<double>& m);

template AVec<float> flatten(const AMat<float>& m);
template AVec<double> flatten(const AMat<double>& m);
template AVec<float> fortran_flatten(const AMat<float>& m);
template AVec<double> fortran_flatten(const AMat<double>& m);
template void diag(const AMat<float>& m, AVec<float> v);
template void diag(const AMat<double>& m, AVec<double> v);



/*---------------------- Linear algebra tools ------------------*/
/*
 * Return x'*M*x
 */
template <class T>
T xtMx(const AVec<T>& x, const AMat<T>& M)
{
    AVec<T> Mx = M*x;
    return x.dot(Mx);
}

template float xtMx(const AVec<float>& x, const AMat<float>& M);
template double xtMx(const AVec<double>& x, const AMat<double>& M);

template <class T>
AMat<T> MtM(const AMat<T>& M)
{
    //
    // BUGBUG, write a faster implementation of this
    // function without the intermediate trasposed() matrix
    //
    AMat<T> Mt = M.transposed();
    return Mt*M;
}

template AMat<float> MtM(const AMat<float>& M);
template AMat<double> MtM(const AMat<double>& M);

/*
 * Use Lapack SVD to invert a matrix
 * Works as psudo-inverse for sigular matrices
 */
AMatd SVDInvert(const AMatd& A)
{
    const double epsilon=1e-12;


    assert(A.rows() == A.columns());

    int d = A.rows();
    AMatd U(d,d);
    AMatd Vt(d,d);
    AVecd vW(d);

    SVD(A, U, vW, Vt);

    AMatd mW(d,d);
    mW.assign_diag(vW);

    U.transpose();
    Vt.transpose();

    // Invert W (the matrix that has values only on the diagonal)
    // for stability: set values on the diagonal to 0 when they are
    // smaller than epsilon.
    for (int i = 0; i < d; ++i)
    {
        if (mW(i,i) < epsilon)
        {
            mW(i,i) = 0;
        }
        else
        {
            mW(i,i) = 1/mW(i,i);
        }
    }

    return Vt*(mW*U);
}

void SVD (
    const AMatd& A,     // matrix to solve
    AMatd& U,           // result U
    AVecd& W,           // result diagonal "matrix" W
    AMatd& Vt           // result V^t
    )
{


//    assert (A.rows() == A.columns());
//    int d = A.rows();
	long m = A.rows();
	long n = A.columns();

    U.resize(m,m); 
    Vt.resize(n,n);


#ifndef USENR

    AVecd fA(fortran_flatten(A));
    AVecd fU(m*m);
    AVecd fVt(n*n);

    long info;
    char job='A';
//    long ld=d;

//    double dwork_size;
    long work_size = -1;

//    std::cout << A << std::endl;

#if 0
    // call dgesvd to compute optimal "work size"
    dgesvd_(
        &job, &job,
        &ld, &ld,
        fA, &ld,
        W,
        fU, &ld,
        fVt, &ld,
        &dwork_size, &work_size,
        &info);

    work_size = long(dwork_size);
    AVecd work(work_size);
#else
    work_size = m*n*20;
    AVecd work(work_size);
#endif

#ifdef USEMKL
#error aoeu
    dgesvd(
        &job, &job,
        (int*)&ld, (int*)&ld,
        fA, (int*)&ld,
        W,
        fU, (int*)&ld,
        fVt, (int*)&ld,
        work, (int*)&work_size,
        (int*)&info);
#else
	long mindim = min(n,m);
    dgesvd_(
        &job, &job,
        &m, &n,
        fA, &m,
        W,
        fU, &m,
        fVt, &n,
        work, &work_size,
        &info);
#endif

    U.assign(fU, true);
    Vt.assign(fVt, true);

#else

    // numerical recipes version


    // copy A into u
    U = A;

    std::vector<double*> nrU(m);
    std::vector<double*> nrVt(n);

    for (int i=0; i<m; i++) {
	nrU[i] = &U(i,0)-1;
    }

    for (int i=0; i<n; i++) {
	nrVt[i] = &Vt(i,0)-1;
    }

    NR::svdcmp<double>(&nrU[0]-1, m, n, &W[0]-1, &nrVt[0]-1);

    Vt.transpose();



#endif
}

/*
 * Solve Ax=b using the above SVD
 */
void SVDSolve(
    const AMatd& A,     // matrix to solve
    const AVecd& b,
    AVecd& x
    )
{
//    assert(A.rows() == A.columns());

//    std::cout << A << std::endl<< std::endl<< std::endl;

    const double epsilon=1e-12;
//    int d = A.rows();
	int m = A.rows();
	int n = A.columns();
    AMatd U(m,m);
    AMatd Vt(n,n);
    AVecd vW(n);

    SVD(A, U, vW, Vt);

    AMatd mW(n,m);
    mW.assign_diag(vW);

    U.transpose();
    Vt.transpose();

    // Invert W (the matrix that has values only on the diagonal)
    // for stability: set values on the diagonal to 0 when they are
    // smaller than epsilon.
	int diaglen = min(n,m);
    for (int i = 0; i < diaglen; ++i)
    {
        if (mW(i,i) < epsilon)
        {
            mW(i,i) = 0;
        }
        else
        {
            mW(i,i) = 1/mW(i,i);
        }
    }

//    std::cout << U << std::endl << vW << std::endl<< std::endl << Vt << std::endl;

    x = Vt*(mW*(U*b));

//    std::cout << x << std::endl;
}

void EigenUpperTriangular(const AMatd& A, AMatd& V, AVecd& evalues)
{
    assert(A.rows() == A.columns());
    long dim = A.rows();
    V.resize(dim, dim);
    evalues.resize(dim);

#ifndef USENR
    AVecd fA(fortran_flatten(A));
    AVecd Work(dim * dim);

    char job = 'V';
    char upl = 'U';
    long work_size = dim*dim;
    long info;
#ifdef USEMKL
    dsyev(&job, &upl, (int*)&dim, fA, (int*)&dim, evalues, Work, (int*)&work_size, (int*)&info);
#else
    dsyev_(&job, &upl, &dim, fA, &dim, evalues, Work, &work_size, &info);
#endif
    assert(info == 0);

    fortran_expand(fA, V);


#else

    // numerical recipes version
    std::cerr<<"EigenUpperTriangular not implemented yet - look at amatc.cpp"<<std::endl;
    exit(-1);

#endif

}

GTB_END_NAMESPACE
