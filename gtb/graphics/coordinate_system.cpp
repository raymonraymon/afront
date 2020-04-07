
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
#include <gtb/graphics/coordinate_system.hpp>
#include <gtb/graphics/line3.hpp>
#include <gtb/graphics/color_rgb.hpp>
#include <gtb/graphics/ogltools.h>
#endif // WIN32


#ifdef OUTLINE
#define inline
#include <gtb/graphics/coordinate_system.ipp>
#undef inline
#endif


GTB_BEGIN_NAMESPACE


template <class T>
tCoordinateSystem<T>::tCoordinateSystem() : 
    m_origin(tPoint3<T>::ZERO),
    m_x(tVector3<T>::VECTOR3_POSITIVE_X),
    m_y(tVector3<T>::VECTOR3_POSITIVE_Y),
    m_z(tVector3<T>::VECTOR3_POSITIVE_Z)
{
}


template <class T>
tCoordinateSystem<T>::tCoordinateSystem(const tCoordinateSystem &cs)
	: m_origin(cs.m_origin),
	  m_x(cs.m_x),
	  m_y(cs.m_y),
	  m_z(cs.m_z)
{
}


template <class T>
tCoordinateSystem<T>::tCoordinateSystem(const Point3 &a_origin,
				   const Vector3 &a_x,
				   const Vector3 &a_y,
				   const Vector3 &a_z)
    : m_origin(a_origin),
      m_x(a_x),
      m_y(a_y),
      m_z(a_z)
{
}


template <class T>
tCoordinateSystem<T> &tCoordinateSystem<T>::operator=(const tCoordinateSystem &cs)
{
    m_origin = cs.m_origin;
    m_x = cs.m_x;
    m_y = cs.m_y;
    m_z = cs.m_z;
    return *this;
}


template <class T>
tMatrix4<T> tCoordinateSystem<T>::matrix() const
{
    return tMatrix4<T>(m_x.x(), m_y.x(), m_z.x(), m_origin.x(),
                   m_x.y(), m_y.y(), m_z.y(), m_origin.y(),
                   m_x.z(), m_y.z(), m_z.z(), m_origin.z(),
                   0.0, 0.0, 0.0, 1.0);
}


template <class T>
tMatrix4<T> tCoordinateSystem<T>::inverse_matrix() const
{
	tMatrix4<T> ir(m_x.x(), m_x.y(), m_x.z(), 0.0,
		   m_y.x(), m_y.y(), m_y.z(), 0.0,
		   m_z.x(), m_z.y(), m_z.z(), 0.0,
		   0.0, 0.0, 0.0, 1.0);

	tMatrix4<T> it(1.0, 0.0, 0.0, -m_origin.x(),
		   0.0, 1.0, 0.0, -m_origin.y(),
		   0.0, 0.0, 1.0, -m_origin.z(),
		   0.0, 0.0, 0.0, 1.0);

	return ir * it;
}


template <class T>
tCoordinateSystem<T> &tCoordinateSystem<T>::rotate(const Line3 &l, value_type theta)
{
    m_origin.rotate(l, theta);
    m_x.rotate(l.direction(), theta);
    m_y.rotate(l.direction(), theta);
    m_z.rotate(l.direction(), theta);
    return *this;
}


template <class T>
tCoordinateSystem<T> &tCoordinateSystem<T>::translate(const Vector3 &t)
{
    m_origin.translate(t);
    return *this;
}


template <class T>
tCoordinateSystem<T> &tCoordinateSystem<T>::move_to(const Point3 &p)
{
    m_origin = p;
    return *this;
}


template <class T>
void tCoordinateSystem<T>::render(value_type length) const
{
    glPushAttrib(GL_ENABLE_BIT);
    glDisable(GL_LIGHTING);

    const Point3 &o = m_origin;
    Point3 px = o + length * m_x;
    Point3 py = o + length * m_y;
    Point3 pz = o + length * m_z;

    COLOR_RGB_RED.load();
    glBegin(GL_LINES);
    o.load();
    px.load();
    glEnd();

    COLOR_RGB_GREEN.load();
    glBegin(GL_LINES);
    o.load();
    py.load();
    glEnd();

    COLOR_RGB_BLUE.load();
    glBegin(GL_LINES);
    o.load();
    pz.load();
    glEnd();

    glPopAttrib();

//  	fprintf(stderr, "o=%f %f %f\n", o.x(), o.y(), o.z());
//  	fprintf(stderr, "x=%f %f %f\n", px.x(), px.y(), px.z());
//  	fprintf(stderr, "y=%f %f %f\n", py.x(), py.y(), py.z());
//  	fprintf(stderr, "z=%f %f %f\n", pz.x(), pz.y(), pz.z());
}

template class tCoordinateSystem<float>;
template class tCoordinateSystem<double>;

GTB_END_NAMESPACE
