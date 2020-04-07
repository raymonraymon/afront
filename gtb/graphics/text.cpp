
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
#include <gtb/graphics/text.hpp>
#include <gtb/graphics/ogltools.h>
#include <stdio.h>
#endif // WIN32


GTB_BEGIN_NAMESPACE



void draw_string(const char *s, int x, int y)
{
	glPushAttrib(GL_CURRENT_BIT);
	glRasterPos2i(x, y);
	while (*s) {
		//glutBitmapCharacter(GLUT_BITMAP_8_BY_13, *s);
		glutBitmapCharacter(GLUT_BITMAP_9_BY_15, *s);
		s++;
	}
	glPopAttrib();
}



void draw_int(int n, int x, int y)
{
	char buf[20];
	sprintf(buf, "%d", n);
	draw_string(buf, x, y);
}


void draw_text_begin(int l, int r, int b, int t)
{
	glPushAttrib(GL_ENABLE_BIT | GL_VIEWPORT_BIT);
	glDisable(GL_LIGHTING);
	glDisable(GL_DEPTH_TEST);

	glViewport(l, b, r - l, t - b);

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(l, r, b, t);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
}


void draw_text_end()
{
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

	glPopAttrib();
}


GTB_END_NAMESPACE
