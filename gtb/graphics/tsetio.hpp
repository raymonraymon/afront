
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


#ifndef __TSETIO_H
#define __TSETIO_H

//
// Triangle set and Point set input / output routines
//
// Copyright (c) 2001 bla bla
//

#include <gtb/common.hpp>
#include <gtb/graphics/indexed_triangle_set.hpp>
#include <gtb/graphics/surfelset.hpp>

GTB_BEGIN_NAMESPACE

void read_obj(FILE* f, IndexedTriangleSet& triangles);
void read_obj(const char* name, IndexedTriangleSet& triangles);
void write_obj(const char* name, const IndexedTriangleSet& triangles, int addidx=0);

void read_triangleset(const char* name, IndexedTriangleSet& triangles);

GTB_END_NAMESPACE

#endif // __TSETIO_H
