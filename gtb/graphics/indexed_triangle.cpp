
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
#ifndef WIN32
#include <algorithm>
#include <gtb/graphics/indexed_triangle.hpp>
#include <gtb/io/io.hpp>
#endif // WIN32


#ifdef OUTLINE
#define inline
#include <gtb/graphics/indexed_triangle.ipp>
#undef inline
#endif


GTB_BEGIN_NAMESPACE


IndexedTriangle::IndexedTriangle()
{
	_indices[0] = -1;
	_indices[1] = -1;
	_indices[2] = -1;
}


IndexedTriangle::IndexedTriangle(int arg_A,
				 int arg_B,
				 int arg_C)
{
	_indices[0] = arg_A;
	_indices[1] = arg_B;
	_indices[2] = arg_C;
}


void IndexedTriangle::read(FILE *fp)
{
	read_int(&_indices[0], fp);
	read_int(&_indices[1], fp);
	read_int(&_indices[2], fp);
}


void IndexedTriangle::write(FILE *fp) const
{
	write_int(_indices[0], fp);
	write_int(_indices[1], fp);
	write_int(_indices[2], fp);
}


void IndexedTriangle::flip()
{
    std::swap(_indices[0], _indices[1]);
}

GTB_END_NAMESPACE
