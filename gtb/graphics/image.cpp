
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
#include <gtb/graphics/image.hpp>
#include <gtb/graphics/ogltools.h>
#include <assert.h>
#endif // WIN32


#ifdef OUTLINE
#define inline
#include <gtb/graphics/image.ipp>
#undef inline
#endif


GTB_BEGIN_NAMESPACE


Image::Image()
	: _pixels(NULL)
{
}


Image::Image(GLsizei a_width, GLsizei a_height, GLenum a_format, GLenum a_type)
	: _width(a_width),
	  _height(a_height),
	  _format(a_format),
	  _type(a_type)
{
	assert(_width > 0);
	assert(_height > 0);
	assert(is_format_valid());
	assert(is_type_valid());
	unsigned nbytes = _width * _height * bits_per_pixel() / 8;
	_pixels = new GLubyte[nbytes];
}


void Image::resize(GLsizei a_width, GLsizei a_height)
{
	assert(_pixels != NULL);
	if ((_width != a_width) || (_height != a_height)) {
		delete[] _pixels;
		_width = a_width;
		_height = a_height;
		unsigned nbytes = _width * _height * bits_per_pixel() / 8;
		_pixels = new GLubyte[nbytes];
	}
}


Image::~Image()
{
	if (_pixels != NULL) {
		delete[] _pixels;
		_pixels = NULL;
	}
}


void Image::draw() const
{
	assert(is_valid());
	glPixelStorei(GL_UNPACK_ROW_LENGTH, _width);
	glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0);
	glPixelStorei(GL_UNPACK_SKIP_ROWS, 0);
	glDrawPixels(_width, _height, _format, _type, _pixels);
}


/*
 * render the image
 *   x,y : lower-left corner of image
 *   w,h : width, height to draw
 *      i.e. x,y,w,h define the window in the image to draw
 *   sx,sy : lower-left corner on screen
 */
void Image::draw(unsigned x, unsigned y, unsigned w, unsigned h, unsigned sx, unsigned sy) const
{
	assert(is_valid());
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1); 
	glPixelStorei(GL_UNPACK_ROW_LENGTH, _width);
	glPixelStorei(GL_UNPACK_SKIP_PIXELS, x);
	glPixelStorei(GL_UNPACK_SKIP_ROWS, y);
    glRasterPos2i(sx, sy);
	glDrawPixels(w, h, _format, _type, _pixels);
}


bool Image::is_format_valid() const
{
	bool ret;

	switch (_format) {
	case GL_COLOR_INDEX:
	case GL_RGB:
	case GL_RGBA:
	case GL_BGR:
	case GL_BGRA:
	case GL_RED:
	case GL_GREEN:
	case GL_BLUE:
	case GL_ALPHA:
	case GL_LUMINANCE:
	case GL_LUMINANCE_ALPHA:
	case GL_STENCIL_INDEX:
	case GL_DEPTH_COMPONENT:
		ret = true;
		break;
	default:
		ret = false;
		break;
	}
	return ret;
}


bool Image::is_type_valid() const
{
	bool ret;

	switch (_type) {
	case GL_UNSIGNED_BYTE:
	case GL_BYTE:
	case GL_BITMAP:
	case GL_UNSIGNED_SHORT:
	case GL_SHORT:
	case GL_UNSIGNED_INT:
	case GL_INT:
	case GL_FLOAT:
	case GL_UNSIGNED_BYTE_3_3_2:
	case GL_UNSIGNED_BYTE_2_3_3_REV:
	case GL_UNSIGNED_SHORT_5_6_5:
	case GL_UNSIGNED_SHORT_5_6_5_REV:
	case GL_UNSIGNED_SHORT_4_4_4_4:
	case GL_UNSIGNED_SHORT_4_4_4_4_REV:
	case GL_UNSIGNED_SHORT_5_5_5_1:
	case GL_UNSIGNED_SHORT_1_5_5_5_REV:
	case GL_UNSIGNED_INT_8_8_8_8:
	case GL_UNSIGNED_INT_8_8_8_8_REV:
	case GL_UNSIGNED_INT_10_10_10_2:
	case GL_UNSIGNED_INT_2_10_10_10_REV:
		ret = true;
		break;
	default:
		ret = false;
		break;
	}
	return ret;
}


bool Image::is_valid() const
{
	return (is_format_valid()
		&& is_type_valid()
		&& (_pixels != NULL)
		&& (_width > 0)
		&& (_height > 0));
}


unsigned Image::num_components() const
{
	unsigned nc = 0;
	switch (_format) {
	case GL_COLOR_INDEX:
		nc = 1;
		break;
	case GL_RGB:
		nc = 3;
		break;
	case GL_RGBA:
		nc = 4;
		break;
	case GL_BGR:
		nc = 3;
		break;
	case GL_BGRA:
		nc = 4;
		break;
	case GL_RED:
		nc = 1;
		break;
	case GL_GREEN:
		nc = 1;
		break;
	case GL_BLUE:
		nc = 1;
		break;
	case GL_ALPHA:
		nc = 1;
		break;
	case GL_LUMINANCE:
		nc = 1;
		break;
	case GL_LUMINANCE_ALPHA:
		nc = 2;
		break;
	case GL_STENCIL_INDEX:
		nc = 1;
		break;
	case GL_DEPTH_COMPONENT:
		nc = 1;
		break;
	default:
		GTB_ERROR("invalid format");
		break;
	}
	assert(nc != 0);
	return nc;
}


unsigned Image::bits_per_pixel() const
{
	assert(is_format_valid());
	assert(is_type_valid());
	unsigned bpp = 0;
	switch (_type) {
// unpacked types
	case GL_UNSIGNED_BYTE:
		bpp = 8 * num_components();
		break;
	case GL_BYTE:
		bpp = 8 * num_components();
		break;
	case GL_BITMAP:
		// Can we have more than 1 component here?
		bpp = 1 * num_components();
		break;
	case GL_UNSIGNED_SHORT:
		bpp = 16 * num_components();
		break;
	case GL_SHORT:
		bpp = 16 * num_components();
		break;
	case GL_UNSIGNED_INT:
		bpp = 32 * num_components();
		break;
	case GL_INT:
		bpp = 32 * num_components();
		break;
	case GL_FLOAT:
		bpp = 8 * sizeof(GLfloat) * num_components();
		break;
// packed types
	case GL_UNSIGNED_BYTE_3_3_2:
		bpp = 8;
		break;
	case GL_UNSIGNED_BYTE_2_3_3_REV:
		bpp = 8;
		break;
	case GL_UNSIGNED_SHORT_5_6_5:
		bpp = 16;
		break;
	case GL_UNSIGNED_SHORT_5_6_5_REV:
		bpp = 16;
		break;
	case GL_UNSIGNED_SHORT_4_4_4_4:
		bpp = 16;
		break;
	case GL_UNSIGNED_SHORT_4_4_4_4_REV:
		bpp = 16;
		break;
	case GL_UNSIGNED_SHORT_5_5_5_1:
		bpp = 16;
		break;
	case GL_UNSIGNED_SHORT_1_5_5_5_REV:
		bpp = 16;
		break;
	case GL_UNSIGNED_INT_8_8_8_8:
		bpp = 32;
		break;
	case GL_UNSIGNED_INT_8_8_8_8_REV:
		bpp = 32;
		break;
	case GL_UNSIGNED_INT_10_10_10_2:
		bpp = 32;
		break;
	case GL_UNSIGNED_INT_2_10_10_10_REV:
		bpp = 32;
		break;
	}
	assert(bpp != 0);
	return bpp;
}


void Image::read_from_frame_buffer(GLint x, GLint y)
{
	glReadPixels(x, y, _width, _height, _format, _type, _pixels);
}


void Image::draw_to_frame_buffer()
{
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(0, _width, 0, _height);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	glRasterPos2i(0, 0);
	glDrawPixels(_width, _height, _format, _type, _pixels);

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
}


void Image::clear()
{
	unsigned nbytes = _width * _height * bits_per_pixel() / 8;
    memset(_pixels, 0, nbytes);
}

GTB_END_NAMESPACE
