
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



#ifndef LLS_WRAPPER_H
#define LLS_WRAPPER_H


#include "common.h"
#include "lsqr.h"

template <typename real_type>
class LLSWrapper {
public:
	LLSWrapper(int nrows, int ncols);
	~LLSWrapper();

	void InsertA(int r, int c, real_type val);
	void InsertB(int r, real_type val);
	void InsertX(int r, real_type val);
	real_type GetX(int r);

	void Solve();

private:

	class IndexValPair {
	public:
		IndexValPair() { }
		IndexValPair(int _i, real_type _v) : i(_i), v(_v) { }
		int i;
		real_type v;
	};

	// If MODE = 0, compute  y = y + A*x,
	void multA(dvec *x, dvec *y);

	// If MODE = 1, compute  x = x + A^T*y.
	void multAt(dvec *x, dvec *y);


	static void static_mult(long mode, dvec *x, dvec *y, void *ptr);


	vector< vector<IndexValPair> > cols;
	vector< vector<IndexValPair> > rows;

	vector<real_type> x;
	vector<real_type> b;
};

#include "lls_wrapper.inl"

#endif
