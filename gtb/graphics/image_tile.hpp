
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


#ifndef GTB_IMAGE_TILE_INCLUDED
#define GTB_IMAGE_TILE_INCLUDED

#include <gtb/common.hpp>
#include <gtb/graphics/image.hpp>

GTB_BEGIN_NAMESPACE


class ImageTile {
public:
	ImageTile(GLint x,
		  GLint y,
		  GLsizei width,
		  GLsizei height,
		  GLenum format,
		  GLenum type);

	~ImageTile();

	void read();
	const Image &image() const;
	Image &image();

protected:
	GLint m_x;
	GLint m_y;
	Image *m_img;
};


GTB_END_NAMESPACE

#ifndef OUTLINE
#include <gtb/graphics/image_tile.ipp>
#endif

#endif // GTB_IMAGE_TILE_INCLUDED
