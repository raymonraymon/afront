
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


#ifndef GTB_MATERIAL_INCLUDED
#define GTB_MATERIAL_INCLUDED

#include <gtb/common.hpp>
#include <gtb/graphics/color_rgb.hpp>


GTB_BEGIN_NAMESPACE


class Material {
public:
	Material();
	~Material();

	void set_ambient_color(const ColorRgb &c);
	void set_diffuse_color(const ColorRgb &c);
	void set_specular_color(const ColorRgb &c);
	void set_emissive_color(const ColorRgb &c);
	void set_shininess(real_t shininess);
	void set_transparency(real_t transparency);

	const ColorRgb &ambient_color() const;
	const ColorRgb &diffuse_color() const;
	const ColorRgb &specular_color() const;
	const ColorRgb &emissive_color() const;
	real_t shininess() const;
	real_t transparency() const;

	void load() const;

protected:
	ColorRgb _ambient_color;
	ColorRgb _diffuse_color;
	ColorRgb _specular_color;
	ColorRgb _emissive_color;
	real_t _shininess;
	real_t _transparency;
};


GTB_END_NAMESPACE

#ifndef OUTLINE
#include <gtb/graphics/material.ipp>
#endif

#endif // GTB_MATERIAL_INCLUDED
