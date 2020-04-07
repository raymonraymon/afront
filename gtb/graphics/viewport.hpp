
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


#ifndef GTB_VIEWPORT_INCLUDED
#define GTB_VIEWPORT_INCLUDED

#include <gtb/common.hpp>
#include <gtb/graphics/point2.hpp>

GTB_BEGIN_NAMESPACE


class Viewport {
public:
	Viewport();
	Viewport(int x_min, int y_min, int width, int height);
	Viewport(const Viewport &vp);
	Viewport &operator=(const Viewport &vp);

	int x_min() const;
	int y_min() const;

	int x_max() const;
	int y_max() const;

	int x_center() const;
	int y_center() const;
	Point2 center() const;

	int width() const;
	int height() const;

	void resize(int x_min, int y_min, int width, int height);

	void load() const;

protected:
	int _x_min, _y_min;
	int _width, _height;
};


GTB_END_NAMESPACE

#ifndef OUTLINE
#include <gtb/graphics/viewport.ipp>
#endif

#endif // GTB_VIEWPORT_INCLUDED
