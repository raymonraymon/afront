
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
#include <gtb/graphics/rtpi.hpp>
#include <gtb/io/io.hpp>
#include <gtb/math/math.hpp>
#endif // WIN32


#ifdef OUTLINE
#define inline
#include <gtb/graphics/rtpi.ipp>
#undef inline
#endif


GTB_BEGIN_NAMESPACE


Rtpi::Rtpi()
	: _r(0.0),
	  _t(0.0),
	  _p(0.0),
	  _i(0)
{
}


Rtpi::Rtpi(real_t arg_r, real_t arg_t, real_t arg_p, int arg_i)
	: _r(arg_r),
	  _t(arg_t),
	  _p(arg_p),
	  _i(arg_i)
{
}


Point3 Rtpi::point() const
{
//  	real_t x = _r * cos(_p) * cos(_t);
//  	real_t y = _r * cos(_p) * sin(_t);
//  	real_t z = _r * sin(_p);

	// x points to the right, y points up, and the scanner is
	// looking down the negative z axis.  Theta is the rotation
	// around y, and phi is the rotation around x.
  	real_t x = -_r * cos(_p) * sin(_t);
  	real_t y = _r * sin(_p);
  	real_t z = -_r * cos(_p) * cos(_t);

	return Point3(x, y, z);
}


GTB_END_NAMESPACE
