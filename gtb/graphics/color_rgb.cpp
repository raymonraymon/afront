
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
#include <gtb/graphics/color_rgb.hpp>
#include <gtb/io/io.hpp>
#endif // WIN32


#ifdef OUTLINE
#define inline
#include <gtb/graphics/color_rgb.ipp>
#undef inline
#endif


GTB_BEGIN_NAMESPACE


const ColorRgb COLOR_RGB_BLACK(0.0, 0.0, 0.0);
const ColorRgb COLOR_RGB_RED(1.0, 0.0, 0.0);
const ColorRgb COLOR_RGB_GREEN(0.0, 1.0, 0.0);
const ColorRgb COLOR_RGB_BLUE(0.0, 0.0, 1.0);
const ColorRgb COLOR_RGB_YELLOW(1.0, 1.0, 0.0);
const ColorRgb COLOR_RGB_CYAN(0.0, 1.0, 1.0);
const ColorRgb COLOR_RGB_MAGENTA(1.0, 0.0, 1.0);
const ColorRgb COLOR_RGB_WHITE(1.0, 1.0, 1.0);
const ColorRgb COLOR_RGB_GRAY25(0.25, 0.25, 0.25);
const ColorRgb COLOR_RGB_GRAY50(0.5, 0.5, 0.5);
const ColorRgb COLOR_RGB_GRAY75(0.75, 0.75, 0.75);


GTB_END_NAMESPACE
