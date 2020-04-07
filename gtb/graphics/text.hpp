
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


#ifndef GTB_TEXT_INCLUDED
#define GTB_TEXT_INCLUDED

#include <gtb/common.hpp>

GTB_BEGIN_NAMESPACE


void draw_string(const char *s, int x, int y);
void draw_int(int n, int x, int y);
void draw_text_begin(int l, int r, int b, int t);
void draw_text_end();


GTB_END_NAMESPACE

#endif // GTB_TEXT_INCLUDED
