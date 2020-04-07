
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


#ifndef GTB_IMAGE_TILE_SET_INCLUDED
#define GTB_IMAGE_TILE_SET_INCLUDED

#include <gtb/common.hpp>
#include <gtb/graphics/image_tile.hpp>
#include <vector>
#include <list>

GTB_BEGIN_NAMESPACE


class ImageTileSet {
public:
	ImageTileSet(GLsizei img_x,
		     GLsizei img_y,
		     GLsizei img_width,
		     GLsizei img_height,
		     GLsizei tile_width,
		     GLsizei tile_height,
		     GLenum format,
		     GLenum type);
	ImageTileSet(const ImageTileSet &tiles);
	~ImageTileSet();

	void resize(GLsizei img_x,
		     GLsizei img_y,
		     GLsizei img_width,
		     GLsizei img_height);

	std::list<ImageTile *> &active_tiles();
	void activate_all();
	void deactivate(std::list<ImageTile *>::iterator pos);

	void read_active();
	void read_all();

protected:
	void alloc_tiles();
	void free_tiles();

	std::vector<ImageTile *> m_tiles;
	std::list<ImageTile *> m_active_tiles;

	GLsizei m_img_x;
	GLsizei m_img_y;
	GLsizei m_img_width;
	GLsizei m_img_height;
	GLsizei m_tile_width;
	GLsizei m_tile_height;
	GLenum m_format;
	GLenum m_type;
};


GTB_END_NAMESPACE

#ifndef OUTLINE
#include <gtb/graphics/image_tile_set.ipp>
#endif

#endif // GTB_IMAGE_TILE_SET_INCLUDED
