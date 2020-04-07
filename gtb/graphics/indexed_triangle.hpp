
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


#ifndef GTB_INDEXED_TRIANGLE_INCLUDED
#define GTB_INDEXED_TRIANGLE_INCLUDED

#include <gtb/common.hpp>
#include <stdio.h>

GTB_BEGIN_NAMESPACE


class IndexedTriangle {
public:
	IndexedTriangle();
	IndexedTriangle(int A, int B, int C);

	int A() const;
	int B() const;
	int C() const;

	void reset(int A, int B, int C);

	bool operator==(const IndexedTriangle &t) const;
	bool operator!=(const IndexedTriangle &t) const;

	int operator[](unsigned i) const;
	int &operator[](unsigned i);

	void read(FILE *fp);
	void write(FILE *fp) const;

    void flip();

    const int* get_indices() const;
    int* get_indices();
protected:
	int _indices[3];
};


GTB_END_NAMESPACE

#ifndef OUTLINE
#include <gtb/graphics/indexed_triangle.ipp>
#endif

#endif // GTB_INDEXED_TRIANGLE_INCLUDED
