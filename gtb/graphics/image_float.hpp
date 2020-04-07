
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


#ifndef GTB_IMAGE_FLOAT_INCLUDED
#define GTB_IMAGE_FLOAT_INCLUDED

//
// "Specialized" image class for real_t type
//
#include <gtb/common.hpp>
#include <gtb/graphics/ogltools.h>
#include <gtb/graphics/color_rgb.hpp>


GTB_BEGIN_NAMESPACE


class ImageFloat : public Image {
public:
    ImageFloat(const Image& image);
    ImageFloat(GLsizei width, GLsizei height);

    ColorRgb& operator() (int x, int y);
    const ColorRgb& operator() (int x, int y) const;

    void add(const ImageFloat& rhs, double weight);
    void add(double value);
    void multiply(double value);
    void copy(const ImageFloat& rhs, int x0, int y0, int w, int h, int x, int y);

protected:
    ImageFloat();
    unsigned _row_size;
};


GTB_END_NAMESPACE

#ifndef OUTLINE
#include <gtb/graphics/image_float.inl>
#endif

#endif // GTB_IMAGE_FLOAT_INCLUDED
