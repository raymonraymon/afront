
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
#include <gtb/math/amatc.h>

//#define USEMKL

#ifndef USENR
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
#else
#include "NR/nr.h"
#endif

#include "amatc.h"


GTB_BEGIN_NAMESPACE

/*------------------ AMatc ---------------------*/

template<class T>
AMatc<T> AMatc<T>::operator *(const AMatc<T>& rhs) const
{
    assert(N == rhs.M);

	assert(N==M);	// this only works for square matrices!!


    AMatc<T> r(M, rhs.N);

    const_pointer cpr = rhs.m; // pointer to the current column in the rhs
    pointer epd = r.m;         // pointer to destination element
    for (int i = 0; i < r.N; ++i)
    {
        const_pointer rpl = m;     // pointer to the current row in the lhs
        for (int j = 0; j < r.M; ++j)
        {
            const_pointer epl = rpl;
            const_pointer epr = cpr;
            *epd = 0;
            for (int k = 0; k < N; ++k)
            {
                *epd += (*epl) * (*epr);
                epl += M;
                ++epr;
            }
            ++epd;
            ++rpl;
        }
        cpr += M;
    }

    return r;
}

template<class T>
AVec<T> AMatc<T>::operator *(const AVec<T>& rhs) const
{
    assert(N == rhs.size());

    AVec<T> r(M);
    pointer pr = m;

    for (int i = 0; i < M; ++i)
    {
        pointer p = pr;
        r(i) = 0;
        for (int j = 0; j < N; ++j)
        {
            r[i] += (*p)*rhs[j];
            p += M;
        }
        ++pr;
    }

    return r;
}

template<class T>
void AMatc<T>::transpose()
{
    if (M==N)
    {
        pointer pr = m+M; // upper-right hand side
        pointer pl = m+1; // lower-left
        for (int i = 0; i < M-1; ++i)
        {
            pointer epr = pr;
            pointer epl = pl;
            for (int j = i+1; j < N; ++j)
            {
                std::swap(*epl, *epr);
                epr += M;
                ++epl;
            }
            pr += M + 1;
            pl += M + 1;
        }
    }
    else
    {
        AMatc temp = transposed();
        *this = temp;
    }
}


template<class T>
AMatc<T> AMatc<T>::transposed() const
{
    int N1 = M;
    int M1 = N;
    AMatc result(M1, N1);
    for (int i = 0; i < M1; ++i)
        for (int j = 0; j < N1; ++j)
            result(i,j) = (*this)(j,i);
    return result;
}

template<class T>
AMatc<T>& AMatc<T>::operator=(const AMatc& rhs)
{
    cleanup();
    M=rhs.M;
    N=rhs.N;
    allocate();
    std::copy(rhs.m, rhs.m+M*N, m);
    return *this;
}

/*
 * If ftn==true, transpose the matrix
 * fotran matrices are flat and transposed
 */
template<class T>
void AMatc<T>::assign(typename AMatc<T>::const_pointer v, bool ftn)
{
    std::copy(v, v+M*N, m);
    if (!ftn) transpose(); // BUGBUG: Performance
}

template<class T>
void AMatc<T>::assign(typename AMatc<T>::const_pointer v, int row, int col)
{
    cleanup();
    M=row;
    N=col;
    allocate();

    assign(v);
}

/*
 * Create a diagonal matrix from a vector
 */
template<class T>
void AMatc<T>::assign_diag(const AVec<T>& v)
{
    pointer p = m;
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < M; ++j)
        {
            if (i == j) *p = v[i];
            else *p = 0;
            ++p;
        }
    }
}

template<class T>
void AMatc<T>::allocate()
{
    m = new value_type[M*N];
}

template<class T>
void AMatc<T>::cleanup()
{
    delete[] m;
}

template<class T>
void AMatc<T>::resize(int row, int col)
{
    if ((row != M) || (col != N))
    {
        cleanup();
        M=row;
        N=col;
        allocate();
    }
}


/*
 * Return the rth row as a vector
 */
template <class T>
void AMatc<T>::get_row(int r, AVec<T>& vrow)
{
    vrow.resize(columns());
    for (int i = 0; i < columns(); ++i)
    {
        vrow(i) = (*this)(r, i);
    }
}

template <class T>
bool AMatc<T>::is_finite() const
{
    int K = N*M;
    for (int i = 0; i < K; ++i)
    {
        if (!isfinite(m[i])) return false;
    }
    return true;
}

/*
 * Return the cth column as a vector
 */
template <class T>
void AMatc<T>::get_column(int c, AVec<T>& vcol)
{
    vcol.resize(rows());
    for (int i = 0; i < rows(); ++i)
    {
        vcol(i) = (*this)(i, c);
    }
}

/*------------------ Utilities -----------------------*/
template <class T>
std::ostream& operator << (std::ostream& s, const AMatc<T>& m)
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
 * Copy the diagonal of a matrix to a vector, assume the size are correct
 */
template <class T>
void diag(const AMatc<T>& m, AVec<T> v)
{
    assert(m.rows() == m.columns() == v.size());
    int N = v.size();
    for (int i = 0; i < N; ++i)
    {
        v(i) = m(i,i);
    }
}

// Instantiate for float and double
template class AMatc<float>;
template class AMatc<double>;
/*
using namespace std;
template ostream& operator << <float> (ostream& s, const AVec<float>& v);
template ostream& operator << <double> (ostream& s, const AVec<double>& v);
template ostream& operator << <float> (ostream& s, const AMatc<float>& m);
template ostream& operator << <double> (ostream& s, const AMatc<double>& m);
*/
template void diag(const AMatc<float>& m, AVec<float> v);
template void diag(const AMatc<double>& m, AVec<double> v);



/*---------------------- Linear algebra tools ------------------*/
/*
 * Return x'*M*x
 */
template <class T>
T xtMx(const AVec<T>& x, const AMatc<T>& M)
{
    AVec<T> Mx = M*x;
    return x.dot(Mx);
}

template float xtMx(const AVec<float>& x, const AMatc<float>& M);
template double xtMx(const AVec<double>& x, const AMatc<double>& M);

template <class T>
AMatc<T> MtM(const AMatc<T>& mat)
{
    AMatc<T> r(mat.N, mat.N);

    typename AMatc<T>::const_pointer cpr = mat.m;  // pointer to the current column in the rhs
    typename AMatc<T>::pointer epd = r.m;        // pointer to destination element
    for (int i = 0; i < r.N; ++i)
    {
        typename AMatc<T>::const_pointer rpl = mat.m;  // pointer to the current row in the lhs
        for (int j = 0; j < r.M; ++j)
        {
            typename AMatc<T>::const_pointer epl = rpl;
            typename AMatc<T>::const_pointer epr = cpr;
            *epd = 0;
            for (int k = 0; k < mat.M; ++k)
            {
                *epd += *epl * (*epr);
                ++epl;
                ++epr;
            }
            ++epd;
            rpl+=mat.M;
        }
        cpr += mat.M;
    }

    return r;
}

template AMatc<float> MtM(const AMatc<float>& M);
template AMatc<double> MtM(const AMatc<double>& M);

/*
 * Use Lapack SVD to invert a matrix
 * Works as psudo-inverse for sigular matrices
 */
template <class T>
AMatc<T> SVDInvert(const AMatc<T>& A)
{
    const double epsilon=1e-12;

    assert(A.rows() == A.columns());

    int d = A.rows();
    AMatc<T> U(d,d);
    AMatc<T> Vt(d,d);
    AVec<T> vW(d);

    SVD(A, U, vW, Vt);

    AMatc<T> mW(d,d);
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

template AMatc<float> SVDInvert(const AMatc<float>& A);
template AMatc<double> SVDInvert(const AMatc<double>& A);


template<>
void SVD<double> (
    const AMatc<double>& A,     // matrix to solve
    AMatc<double>& U,           // result U
    AVec<double>& W,           // result diagonal "matrix" W
    AMatc<double>& Vt           // result V^t
    )
{
    assert (A.rows() == A.columns());
    int d = A.rows();

//    AVec<T> fA(fortran_flatten(A));
//    AVec<T> fU(d*d);
//    AVec<T> fVt(d*d);
    long info;
    char job='A';
    long ld=d;

    U.resize(d,d);
    Vt.resize(d,d);

#ifndef USENR
//    double dwork_size;
    long work_size = -1;

//    std::cout << A << std::endl;

#if 0
    // call dgesvd to compute optimal "work size"
    dgesvd_(
        &job, &job,
        &ld, &ld,
        A.flatten(), &ld,
        W,
        fU, &ld,
        fVt, &ld,
        &dwork_size, &work_size,
        &info);

    work_size = long(dwork_size);
    AVec<double> work(work_size);
#else
    work_size = d*d*8; 
    AVec<double> work(work_size);
#endif

    AVec<double> awork(A.rows() * A.columns(), (double*)A.flatten());

#ifdef USEMKL
    dgesvd( // NOTE dgesdd is slower for our case!!
        &job, &job,
        (int*)&ld, (int*)&ld,
        awork, (int*)&ld,
        W,
        U.flatten(), (int*)&ld,
        Vt.flatten(), (int*)&ld,
        work, (int*)&work_size,
        (int*)&info);
#else
        dgesvd_(
        &job, &job,
        &ld, &ld,
        awork, &ld,
        W,
		U.flatten(), &ld,
		Vt.flatten(), &ld,
        work, &work_size,
        &info);
#endif

#else

    // transpose A into u
    for (int r=0; r<d; r++) {
	for (int c=0; c<d; c++) {
	    U(r,c) = A(c,r);
	}
    }

    std::vector<double*> nrU(d);
    std::vector<double*> nrVt(d);

    for (int i=0; i<d; i++) {
	nrU[i] = &U(0,i)-1;
	nrVt[i] = &Vt(0,i)-1;
    }

    NR::svdcmp<double>(&nrU[0]-1, d, d, &W[0]-1, &nrVt[0]-1);

    U.transpose();

#endif

//    U.resize(d,d); U.assign(fU, true);
//    Vt.resize(d,d); Vt.assign(fVt, true);
}

template<>
void SVD<float> (
    const AMatc<float>& A,     // matrix to solve
    AMatc<float>& U,           // result U
    AVec<float>& W,           // result diagonal "matrix" W
    AMatc<float>& Vt           // result V^t
    )
{

    assert (A.rows() == A.columns());
    int d = A.rows();

//    AVec<T> fA(fortran_flatten(A));
//    AVec<T> fU(d*d);
//    AVec<T> fVt(d*d);
    long info;
    char job='A';
    long ld=d;

    U.resize(d,d);
    Vt.resize(d,d);


#ifndef USENR

//    float dwork_size;
    long work_size = -1;

//    std::cout << A << std::endl;

#if 0
    // call dgesvd to compute optimal "work size"
    sgesvd_(
        &job, &job,
        &ld, &ld,
        A.flatten(), &ld,
        W,
        fU, &ld,
        fVt, &ld,
        &dwork_size, &work_size,
        &info);

    work_size = long(dwork_size);
    AVec<float> work(work_size);
#else
    work_size = d*d*8; 
    AVec<float> work(work_size);
#endif

    AVec<float> awork(A.rows() * A.columns(), (float*)A.flatten());

#ifdef USEMKL
    sgesvd( // NOTE dgesdd is slower for our case!!
        &job, &job,
        (int*)&ld, (int*)&ld,
        awork, (int*)&ld,
        W,
        U.flatten(), (int*)&ld,
        Vt.flatten(), (int*)&ld,
        work, (int*)&work_size,
        (int*)&info);
#else
    sgesvd_(
        &job, &job,
        &ld, &ld,
        awork, &ld,
        W,
        U.flatten(), &ld,
        Vt.flatten(), &ld,
        work, &work_size,
        &info);
#endif


#else

    // transpose A into u
    for (int r=0; r<d; r++) {
	for (int c=0; c<d; c++) {
	    U(r,c) = A(c,r);
	}
    }

    std::vector<float*> nrU(d);
    std::vector<float*> nrVt(d);

    for (int i=0; i<d; i++) {
	nrU[i] = &U(0,i)-1;
	nrVt[i] = &Vt(0,i)-1;
    }

    NR::svdcmp<float>(&nrU[0]-1, d, d, &W[0]-1, &nrVt[0]-1);

    U.transpose();

#endif

//    U.resize(d,d); U.assign(fU, true);
//    Vt.resize(d,d); Vt.assign(fVt, true);
}


/*
 * Solve Ax=b using the above SVD
 */
template <class T>
void SVDSolve(
    const AMatc<T>& A,     // matrix to solve
    const AVec<T>& b,
    AVec<T>& x
    )
{
    assert(A.rows() == A.columns());

//    std::cout << A << std::endl<< std::endl<< std::endl;

    const T epsilon=sizeof(T) == 4 ? 1e-8 : 1e-12; // BUGBUG HACK

    int d = A.rows();
    AMatc<T> U(d,d);
    AMatc<T> Vt(d,d);
    AVec<T> vW(d);

    SVD(A, U, vW, Vt);

    AMatc<T> mW(d,d);
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

    x = Vt*(mW*(U*b));
}

template void SVDSolve<float>(const AMatc<float>& A, const AVec<float>& b, AVec<float>& x);
template void SVDSolve<double>(const AMatc<double>& A, const AVec<double>& b, AVec<double>& x);


void EigenUpperTriangular(AMatc<double>& A, AVec<double>& evalues)
{
    assert(A.rows() == A.columns());
    long dim = A.rows();
    evalues.resize(dim);

#ifndef USENR

//    AVec<double> fA(fortran_flatten(A));
    long work_size = dim*(dim+2);
    AVec<double> Work(work_size);

    char job = 'V';
    char upl = 'U';
    long info;
//    dsyev_(&job, &upl, &dim, fA, &dim, evalues, Work, &work_size, &info);
#ifdef USEMKL
    dsyev(&job, &upl, (int*)&dim, A.flatten(), (int*)&dim, evalues, Work, (int*)&work_size, (int*)&info);
#else
    dsyev_(&job, &upl, &dim, A.flatten(), &dim, evalues, Work, &work_size, &info);
#endif

    assert(info == 0);

#else

    // make sure the whole matrix is filled in
    for (int r=1; r<dim; r++) {
	for (int c=0; c<r; c++) {
	    A(r,c) = A(c,r);
	}
    }

    AMatc<double> nV(dim,dim);

    std::vector<double*> nrA(dim);
    std::vector<double*> nrEvec(dim);
    for (int i=0; i<dim; i++) {
	nrA[i] = &A(0,i)-1;
	nrEvec[i] = &nV(0,i)-1;
    }

    int nrot=0;
    NR::jacobi<double>(&nrA[0]-1, dim, &evalues[0]-1, &nrEvec[0]-1, &nrot);


    // transpose the eigenvectors out of nV
    for (int r=0; r<dim; r++) {
	for (int c=0; c<dim; c++) {
	    A(r,c) = nV(c,r);
	}
    }


    // sort them into increasing order
    for (int first=0; first<dim-1; first++) {
	int min=first;
	for (int i=first+1; i<dim; i++) {
	    if (evalues[i]<evalues[min])
		min = i;
	}

	if (min != first) {
	    std::swap(evalues[first], evalues[min]);

	    for (int i=0; i<dim; i++) {
		std::swap(A(i,min), A(i,first));
	    }
	}
    }

#endif

}



void EigenUpperTriangular(AMatc<float>& A, AVec<float>& evalues)
{
    assert(A.rows() == A.columns());
    long dim = A.rows();
    evalues.resize(dim);

#ifndef USENR

//    AVec<float> fA(fortran_flatten(A));
    long work_size = dim*(dim+2);
    AVec<float> Work(work_size);

    char job = 'V';
    char upl = 'U';
    long info;
//    ssyev_(&job, &upl, &dim, fA, &dim, evalues, Work, &work_size, &info);
#ifdef USEMKL
    ssyev(&job, &upl, (int*)&dim, A.flatten(), (int*)&dim, evalues, Work, (int*)&work_size, (int*)&info);
#else
    ssyev_(&job, &upl, &dim, A.flatten(), &dim, evalues, Work, &work_size, &info);
#endif

    assert(info == 0);

#else

    // make sure the whole matrix is filled in
    for (int r=1; r<dim; r++) {
	for (int c=0; c<r; c++) {
	    A(r,c) = A(c,r);
	}
    }

    AMatc<float> nV(dim,dim);

    std::vector<float*> nrA(dim);
    std::vector<float*> nrEvec(dim);
    for (int i=0; i<dim; i++) {
	nrA[i] = &A(0,i)-1;
	nrEvec[i] = &nV(0,i)-1;
    }

    int nrot=0;
    NR::jacobi<float>(&nrA[0]-1, dim, &evalues[0]-1, &nrEvec[0]-1, &nrot);


    // transpose the eigenvectors out of nV
    for (int r=0; r<dim; r++) {
	for (int c=0; c<dim; c++) {
	    A(r,c) = nV(c,r);
	}
    }


    // sort them into increasing order
    for (int first=0; first<dim-1; first++) {
	int min=first;
	for (int i=first+1; i<dim; i++) {
	    if (evalues[i]<evalues[min])
		min = i;
	}

	if (min != first) {
	    std::swap(evalues[first], evalues[min]);

	    for (int i=0; i<dim; i++) {
		std::swap(A(i,min), A(i,first));
	    }
	}
    }

#endif

}

GTB_END_NAMESPACE
