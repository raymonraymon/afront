
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



#ifndef __GTB_LINEARSOLVERS_H
#define __GTB_LINEARSOLVERS_H


#include <gtb/common.hpp>

GTB_BEGIN_NAMESPACE


// check check relative error to x=[1,1,....,1]
template <typename T, int N>
T relerror(T x[N]) {
	T relerror=0;
	for (int i=0; i<N; i++) {
		relerror += (x[i]-1)*(x[i]-1);
	}
	relerror = sqrt(relerror);
	return (relerror/sqrt((T)N));
}

// solve an upper triangular system
template <typename T, int N>
void bsa(const T A[N][N], const T b[N], T x[N]) {
	for (int i=0; i<N; i++) {
		x[i] = b[i];
	}

	for (int i=N-1; i>=0; i--) {
		for (int j=0; j<i; j++) {
			x[j] -= x[i] * A[j][i] / A[i][i];
		}
		x[i] /= A[i][i];
	}
}

// solve a lower triangular system
template <typename T, int N>
void fsa(const T A[N][N], const T b[N], T x[N]) {
	for (int i=0; i<N; i++) {
		x[i] = b[i];
	}

	for (int i=0; i<N; i++) {
		for (int j=i+1; j<N; j++) {
			x[j] -= x[i] * A[j][i] / A[i][i];
		}
		x[i] /= A[i][i];
	}
}

template <typename T, int N>
void cholesky(const T A[N][N], const T b[N], T x[N]) {

	T l[N][N];
	for (int i=0; i<N; i++) {
		for (int j=0; j<N; j++) {
			l[i][j]=0;
		}
	}

	// first col of l
	l[0][0] = sqrt(A[0][0]);
	for (int i=1; i<N; i++) {
		l[i][0] = A[i][0] / l[0][0];
	}

	// the rest of l
	for (int i=1; i<N; i++) {
		l[i][i] = A[i][i];
		for (int j=0; j<i; j++) {
			l[i][i] -= l[i][j] * l[i][j];
		}
		l[i][i] = sqrt(l[i][i]);

		for (int j=i+1; j<N; j++) {
			l[j][i] = A[j][i];
			for (int k=0; k<i; k++) {
				l[j][i] -= l[j][k]*l[i][k];	
			}
			l[j][i] /= l[i][i];
		}
	}

	// make l transpose
	T lt[N][N];
	for (int i=0; i<N; i++) {
		for (int j=0; j<N; j++) {
			lt[j][i] = l[i][j];
		}
	}

	// solve ly=b
	T y[N];
	fsa<T,N>(l, b, y);

	// solve lt x=y
	bsa<T,N>(lt, y, x);
}

template <typename T, int N>
void householder(const T A[N][N], const T b[N], T x[N]) {

	T r[N][N];
	T w[N][N];	// lower triangular rep. of the H's
	for (int i=0; i<N; i++) {
		for (int j=0; j<N; j++) {
			r[i][j]=A[i][j];
			w[i][j]=0;
		}
	}

	for (int h=0; h<N; h++) {
		T magx = 0;
		for (int i=h; i<N; i++) {
			magx += r[i][h]*r[i][h];
		}
		magx = sqrt(magx);
		T sign = (r[h][h]>0) ? 1.0 : -1.0;
		T k = sqrt(2.0 * magx * (magx + sign*r[h][h]));
		T kinv = 1.0/k;
		for (int i=h; i<N; i++) {
			w[i][h] = r[i][h] * kinv;
		}
		w[h][h] += sign*magx * kinv;

		// apply this w to the lower (N-h)x(N-h) rect of r
		for (int h2=h; h2<N; h2++) {
			T dot = 0;
			for (int i=h; i<N; i++) {
				dot += w[i][h] * r[i][h2];
			}
			for (int i=h; i<N; i++) {
				r[i][h2] -= 2.0*w[i][h]*dot;
			}
		}
	}


	// solve Hy=b
	T y[N];
	for (int i=0; i<N; i++) {
		y[i] = b[i];
	}

	// apply the H's
	for (int h=0; h<N; h++) {
		T dot = 0;
		for (int i=h; i<N; i++) {
			dot += w[i][h] * y[i];
		}
		for (int i=h; i<N; i++) {
			y[i] -= 2.0*w[i][h]*dot;
		}
	}

	// solve Rx=y
	bsa<T,N>(r, y, x);
}

GTB_END_NAMESPACE

#endif


