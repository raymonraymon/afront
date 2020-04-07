
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


#ifndef GTB_LIGHT_INCLUDED
#define GTB_LIGHT_INCLUDED

#include <gtb/graphics/ogltools.h>
#include <gtb/common.hpp>

GTB_BEGIN_NAMESPACE


class Light {
public:
	Light(GLenum light);

	void enable();
	void disable();
	bool is_enabled();

	void set_ambient(GLfloat r, GLfloat g, GLfloat b, GLfloat a);
	void set_diffuse(GLfloat r, GLfloat g, GLfloat b, GLfloat a);
	void set_specular(GLfloat r, GLfloat g, GLfloat b, GLfloat a);
	void set_position(GLfloat x, GLfloat y, GLfloat z, GLfloat w);
	void load() const;

protected:
	GLenum _light;
	bool _is_enabled;
	GLfloat _ambient[4];
	GLfloat _diffuse[4];
	GLfloat _specular[4];
	GLfloat _position[4];
};


GTB_END_NAMESPACE

#ifndef OUTLINE
#include <gtb/graphics/light.ipp>
#endif

#endif // GTB_LIGHT_INCLUDED
