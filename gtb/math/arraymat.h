
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


#ifndef __ARRAYMAT_H
#define __ARRAYMAT_H

GTB_BEGIN_NAMESPACE

/*
 * A wrapper that holds a pointer to a 1d array and treats it as a 2d array
 * Similar to image() with no ownership
 *
 * NOTE: this class is not the owner of the pointer
 */
template <class T>
class ArrayMat
{
public:
	ArrayMat(T* a, unsigned w, unsigned h) :
	  _a(a), W(w), H(h)
	{
	}

	T& operator()(unsigned y, unsigned x)
	{
		assert(y<H);
		assert(x<W);
		return _a[y*W+x];
	}

	T& operator()(unsigned idx)
	{
		assert(idx<W*H);
		return _a[idx];
	}

protected:
	T* _a;
	unsigned W, H;
};

GTB_END_NAMESPACE

#endif // __ARRAYMAT_H
