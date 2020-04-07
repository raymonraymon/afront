
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
#include <gtb/graphics/sphere.hpp>
#endif // WIN32


#ifdef OUTLINE
#define inline
#include <gtb/graphics/sphere.ipp>
#undef inline
#endif


GTB_BEGIN_NAMESPACE


template <class T>
tSphere<T>::tSphere()
{
}


template <class T>
tSphere<T>::tSphere(const Point3 &c, value_type r)
	: _center(c),
	  _radius(r)
{
}


template <class T>
tSphere<T>::tSphere(const tSphere &s)
	: _center(s._center),
	  _radius(s._radius)
{
}


template <class T>
void tSphere<T>::render(int slices,
		    int stacks,
		    GLenum draw_style,
		    GLenum normal_type) const
{
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	const Point3 &so = center();
	glTranslatef(so.x(), so.y(), so.z());
	GLUquadric *qobj = gluNewQuadric();
	gluQuadricDrawStyle(qobj, draw_style);
	gluQuadricNormals(qobj, normal_type);
	gluSphere(qobj, radius(), slices, stacks);
	gluDeleteQuadric(qobj);
	glPopMatrix();
}

template class tSphere<float>;
template class tSphere<double>;

GTB_END_NAMESPACE
