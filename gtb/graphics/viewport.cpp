
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
#include <gtb/graphics/viewport.hpp>
#include <gtb/graphics/plane.hpp>
#include <gtb/graphics/polygon3.hpp>
#include <gtb/graphics/ogltools.h>
#endif // WIN32


#ifdef OUTLINE
#define inline
#include <gtb/graphics/viewport.ipp>
#undef inline
#endif


GTB_BEGIN_NAMESPACE


Viewport::Viewport()
	: _x_min(0),
	  _y_min(0),
	  _width(512),
	  _height(512)
{
}


Viewport::Viewport(int xmin, int ymin, int w, int h)
	: _x_min(xmin),
	  _y_min(ymin),
	  _width(w),
	  _height(h)
{
}


Viewport::Viewport(const Viewport &vp)
	: _x_min(vp._x_min),
	  _y_min(vp._y_min),
	  _width(vp._width),
	  _height(vp._height)
{
}


GTB_END_NAMESPACE
