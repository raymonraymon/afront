
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
#include <gtb/graphics/image_tile_set.hpp>
#endif // WIN32

#ifdef OUTLINE
#define inline
#include <gtb/graphics/image_tile_set.ipp>
#undef inline
#endif

GTB_BEGIN_NAMESPACE


ImageTileSet::ImageTileSet(GLsizei a_img_x,
			   GLsizei a_img_y,
			   GLsizei a_img_width,
			   GLsizei a_img_height,
			   GLsizei a_tile_width,
			   GLsizei a_tile_height,
			   GLenum a_format,
			   GLenum a_type)
	: m_img_x(a_img_x),
	  m_img_y(a_img_y),
	  m_img_width(a_img_width),
	  m_img_height(a_img_height),
	  m_tile_width(a_tile_width),
	  m_tile_height(a_tile_height),
	  m_format(a_format),
	  m_type(a_type)
{
	alloc_tiles();
}


ImageTileSet::ImageTileSet(const ImageTileSet &tiles)
	: m_img_x(tiles.m_img_x),
	  m_img_y(tiles.m_img_y),
	  m_img_width(tiles.m_img_width),
	  m_img_height(tiles.m_img_height),
	  m_tile_width(tiles.m_tile_width),
	  m_tile_height(tiles.m_tile_height),
	  m_format(tiles.m_format),
	  m_type(tiles.m_type)
{
	alloc_tiles();
}


void ImageTileSet::alloc_tiles()
{
	m_tiles.clear();

	// for each row
	for (GLint y = 0; y < m_img_height; y+= m_tile_height) {

		// determine height
		GLsizei h;
		if (y + m_tile_height <= m_img_height) {
			h = m_tile_height;
		} else {
			h = m_img_height - y;
		}
		assert((y + h) <= m_img_height);

		// for each column
		for (GLint x = 0; x < m_img_width; x+= m_tile_width) {

			// determine width
			GLsizei w;
			if (x + m_tile_width <= m_img_width) {
				w = m_tile_width;
			} else {
				w = m_img_width - x;
			}
			assert((x + w) <= m_img_width);

			// create tile
			ImageTile *tile = new ImageTile(x + m_img_x,
							y + m_img_y,
							w, h,
							m_format, m_type);
			m_tiles.push_back(tile);
		}
	}
}


ImageTileSet::~ImageTileSet()
{
	free_tiles();
}


void ImageTileSet::free_tiles()
{
	for (unsigned i = 0; i < m_tiles.size(); i++) {
		delete m_tiles[i];
	}
	m_tiles.clear();
	m_active_tiles.clear();
}


void ImageTileSet::resize(GLsizei a_img_x,
			   GLsizei a_img_y,
			   GLsizei a_img_width,
			   GLsizei a_img_height)
{
	if ((a_img_x != m_img_x)
	    || (a_img_y != m_img_y)
	    || (a_img_width != m_img_width)
	    || (a_img_height != m_img_height)) {
		free_tiles();
		m_img_x = a_img_x;
		m_img_y = a_img_y;
		m_img_width = a_img_width;
		m_img_height = a_img_height;
		alloc_tiles();
	}
}


GTB_END_NAMESPACE
