
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


#include "stdafx.h"
#include "mlslibdefs.h"
#include "poly.h"
#include "poly1.h"
#include "poly2.h"

using namespace gtb;

MLSLIB_BEGIN_NAMESPACE

template <class REAL>
void write_poly(FILE* f, const Poly<REAL>* c)
{
	write_int(c->degree(), f);
	c->write(f);
}

template <class REAL>
Poly<REAL>* read_poly(FILE* f)
{
	int degree;
	read_int(&degree, f);
	Poly<REAL>* c;
	switch (degree)
	{
	case 1:
		c = new Poly1<REAL>;
		break;
	case 2:
		c = new Poly2<REAL>;
		break;
	default:
		assert(0);
	}
	c->read(f);
	return c;
}

template void write_poly(FILE* f, const Poly<float>* c);
template void write_poly(FILE* f, const Poly<double>* c);
template Poly<float>* read_poly(FILE* f);
template Poly<double>* read_poly(FILE* f);

MLSLIB_END_NAMESPACE
