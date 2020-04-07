
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


#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

namespace NR {

// Internal routines
inline void* nralloc(size_t size)
{
//    return _aligned_malloc(size, 8);
    return malloc(size);
}

inline void nrfree(void* p)
{
//    _aligned_free(p);
    free(p);
}


template <class T> 
inline T SQR(const T& x) { return x*x; }


inline double DMAX(double a, double b)
{
	return (a > b) ? a : b;
}

inline double DMIN(double a, double b)
{
	return (a < b) ? a : b;
}


template <typename T>
inline T FMAX(T a, T b)
{
	return (a > b) ? a : b;
}

inline double FMIN(double a, double b)
{
	return (a < b) ? a : b;
}

inline float FMIN(float a, float b)
{
	return (a < b) ? a : b;
}

inline long LMAX(long a, long b)
{
	return (a > b) ? a : b;
}

inline long LMIN(long a, long b)
{
	return (a < b) ? a : b;
}

inline int IMAX(int a, int b)
{
	return (a > b) ? a : b;
}

inline int IMIN(int a, int b)
{
	return (a < b) ? a : b;
}

template <typename T>
inline T SIGN(T a, T b);

template <>
inline double SIGN(double a, double b)
{
	return (b >= 0.0) ? fabs(a) : -fabs(a);
}

template <>
inline float SIGN(float a, float b)
{
	return (b >= 0.0) ? fabsf(a) : -fabsf(a);
}

#if defined(__STDC__) || defined(ANSI) || defined(NRANSI) || defined(__cplusplus) /* ANSI */

void nrerror(const char error_text[]);
float *vector(long nl, long nh);
int *ivector(long nl, long nh);
unsigned char *cvector(long nl, long nh);
unsigned long *lvector(long nl, long nh);
double *dvector(long nl, long nh);
float **matrix(long nrl, long nrh, long ncl, long nch);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);
float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl);
float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch);
float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_vector(float *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_cvector(unsigned char *v, long nl, long nh);
void free_lvector(unsigned long *v, long nl, long nh);
void free_dvector(double *v, long nl, long nh);
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch);
void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch);
void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh);

void print_matrix(float** m, long nrl, long nrh, long ncl, long nch);
void print_vector(float *v, long nl, long nh);

//
// Templated stuff
//
template <class T> T *avector(long nl, long nh);
template <class T> void free_avector(T *v, long nl, long nh);
template <class T> T **dmatrix(long nrl, long nrh, long ncl, long nch);
template <class T> void free_dmatrix(T **m, long nrl, long nrh, long ncl, long nch);

#define NR_END 1
template <class T>
T *avector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
    T *v;
    
    v=(T *)nralloc((size_t) ((nh-nl+1+NR_END)*sizeof(T)));
    if (!v) nrerror("allocation failure in avector()");
    return v-nl+NR_END;
}

template <class T>
void free_avector(T *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
    nrfree((void*) (v+nl-NR_END));
}

template <class T>
T **amatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
    T **m;

    /* allocate pointers to rows */
    m=(T **) nralloc((size_t)((nrow+NR_END)*sizeof(T*)));
    if (!m) nrerror("allocation failure 1 in matrix()");
    m += NR_END;
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl]=(T *) nralloc((size_t)((nrow*ncol+NR_END)*sizeof(T)));
    if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

    /* return pointer to array of pointers to rows */
    return m;
}

template <class T>
void free_amatrix(T **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
    nrfree((void*) (m[nrl]+ncl-NR_END));
    nrfree((void*) (m+nrl-NR_END));
}

#else /* ANSI */
/* traditional - K&R */

void nrerror();
float *vector();
float **matrix();
float **submatrix();
float **convert_matrix();
float ***f3tensor();
double *dvector();
double **dmatrix();
int *ivector();
int **imatrix();
unsigned char *cvector();
unsigned long *lvector();
void free_vector();
void free_dvector();
void free_ivector();
void free_cvector();
void free_lvector();
void free_matrix();
void free_submatrix();
void free_convert_matrix();
void free_dmatrix();
void free_imatrix();
void free_f3tensor();

#endif /* ANSI */
}

#endif /* _NR_UTILS_H_ */
