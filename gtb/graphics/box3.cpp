
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

#include <gtb/graphics/box3.hpp>
#include <gtb/graphics/plane.hpp>
#include <gtb/graphics/polygon3.hpp>
#include <gtb/graphics/ogltools.h>




#ifdef OUTLINE
#define inline
#include <gtb/graphics/box3.ipp>
#undef inline
#endif


using namespace std;


GTB_BEGIN_NAMESPACE


// see AT&T notes, p. 34
template<class T>
const unsigned tBox3<T>::_vertex_indices[6][4] = {
	{0, 1, 3, 2},
	{4, 6, 7, 5},
	{0, 4, 5, 1},
	{2, 3, 7, 6},
	{0, 2, 6, 4},
	{1, 5, 7, 3}
};



template<class T>
typename tBox3<T>::value_type tBox3<T>::shortest_axis_length() const
{
	value_type dx = x_length();
	value_type dy = y_length();
	value_type dz = z_length();
	return min(min(dx, dy), dz);
}


template<class T>
typename tBox3<T>::value_type tBox3<T>::longest_axis_length() const
{
	value_type dx = x_length();
	value_type dy = y_length();
	value_type dz = z_length();
	return max(max(dx, dy), dz);
}


template<class T>
int tBox3<T>::classify_position(const tBox3 &b) const
{
	if ((_max_pt[0] < b._min_pt[0]) || (_min_pt[0] > b._max_pt[0]) ||
	    (_max_pt[1] < b._min_pt[1]) || (_min_pt[1] > b._max_pt[1]) ||
	    (_max_pt[2] < b._min_pt[2]) || (_min_pt[2] > b._max_pt[2])) {
		return OUTSIDE;
	}

	if ((_min_pt[0] <= b._min_pt[0]) && (_max_pt[0] >= b._max_pt[0]) &&
	    (_min_pt[1] <= b._min_pt[1]) && (_max_pt[1] >= b._max_pt[1]) &&
	    (_min_pt[2] <= b._min_pt[2]) && (_max_pt[2] >= b._max_pt[2])) {
		return INSIDE;
        }

	return INTERSECT;
}


template<class T>
tPolygon3<T> tBox3<T>::face(unsigned i) const
{
	assert(i < 6);
	return tPolygon3<T>(vertex(i, 0),
			vertex(i, 1),
			vertex(i, 2),
			vertex(i, 3),
			normal(i));
}


template<class T>
tPlane<T> tBox3<T>::plane(unsigned i) const
{
	assert(i < 6);
	return tPlane<T>(vertex(_vertex_indices[i][0]),
		     vertex(_vertex_indices[i][1]),
		     vertex(_vertex_indices[i][2]));
}


template<class T>
void tBox3<T>::render() const
{
	glBegin(GL_QUADS);
	for (unsigned i = 0; i < 6; i++) {
		normal(i).load_as_normal();
		vertex(i, 0).load();
		vertex(i, 1).load();
		vertex(i, 2).load();
		vertex(i, 3).load();
	}
	glEnd();
}


template<>
void tBox3<float>::outline() const
{
	glPushAttrib(GL_ENABLE_BIT);
	glDisable(GL_LIGHTING);

	glBegin(GL_LINE_LOOP);
	glVertex3f(_min_pt[0], _min_pt[1], _min_pt[2]);
	glVertex3f(_max_pt[0], _min_pt[1], _min_pt[2]);
	glVertex3f(_max_pt[0], _max_pt[1], _min_pt[2]);
	glVertex3f(_min_pt[0], _max_pt[1], _min_pt[2]);
	glEnd();

	glBegin(GL_LINE_LOOP);
	glVertex3f(_min_pt[0], _min_pt[1], _max_pt[2]);
	glVertex3f(_max_pt[0], _min_pt[1], _max_pt[2]);
	glVertex3f(_max_pt[0], _max_pt[1], _max_pt[2]);
	glVertex3f(_min_pt[0], _max_pt[1], _max_pt[2]);
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(_min_pt[0], _min_pt[1], _min_pt[2]);
	glVertex3f(_min_pt[0], _min_pt[1], _max_pt[2]);
	glVertex3f(_max_pt[0], _min_pt[1], _min_pt[2]);
	glVertex3f(_max_pt[0], _min_pt[1], _max_pt[2]);
	glVertex3f(_max_pt[0], _max_pt[1], _min_pt[2]);
	glVertex3f(_max_pt[0], _max_pt[1], _max_pt[2]);
	glVertex3f(_min_pt[0], _max_pt[1], _min_pt[2]);
	glVertex3f(_min_pt[0], _max_pt[1], _max_pt[2]);
	glEnd();

	glPopAttrib();
}

template<>
void tBox3<double>::outline() const
{
	glPushAttrib(GL_ENABLE_BIT);
	glDisable(GL_LIGHTING);

	glBegin(GL_LINE_LOOP);
	glVertex3d(_min_pt[0], _min_pt[1], _min_pt[2]);
	glVertex3d(_max_pt[0], _min_pt[1], _min_pt[2]);
	glVertex3d(_max_pt[0], _max_pt[1], _min_pt[2]);
	glVertex3d(_min_pt[0], _max_pt[1], _min_pt[2]);
	glEnd();

	glBegin(GL_LINE_LOOP);
	glVertex3d(_min_pt[0], _min_pt[1], _max_pt[2]);
	glVertex3d(_max_pt[0], _min_pt[1], _max_pt[2]);
	glVertex3d(_max_pt[0], _max_pt[1], _max_pt[2]);
	glVertex3d(_min_pt[0], _max_pt[1], _max_pt[2]);
	glEnd();

	glBegin(GL_LINES);
	glVertex3d(_min_pt[0], _min_pt[1], _min_pt[2]);
	glVertex3d(_min_pt[0], _min_pt[1], _max_pt[2]);
	glVertex3d(_max_pt[0], _min_pt[1], _min_pt[2]);
	glVertex3d(_max_pt[0], _min_pt[1], _max_pt[2]);
	glVertex3d(_max_pt[0], _max_pt[1], _min_pt[2]);
	glVertex3d(_max_pt[0], _max_pt[1], _max_pt[2]);
	glVertex3d(_min_pt[0], _max_pt[1], _min_pt[2]);
	glVertex3d(_min_pt[0], _max_pt[1], _max_pt[2]);
	glEnd();

	glPopAttrib();
}


template<class T>
tBox3<T> tBox3<T>::bounding_box(const vector<tPoint3<T> > &v)
{
	assert(v.size() > 0);
	tPoint3<T> p = v[0];
	value_type min_x = p.x(), min_y = p.y(), min_z = p.z();
	value_type max_x = p.x(), max_y = p.y(), max_z = p.z();
	for (unsigned i = 1; i < v.size(); i++) {
		p = v[i];
		min_x = min(min_x, p.x());
		min_y = min(min_y, p.y());
		min_z = min(min_z, p.z());
		max_x = max(max_x, p.x());
		max_y = max(max_y, p.y());
		max_z = max(max_z, p.z());
	}
	return tBox3(min_x, min_y, min_z,
		    max_x, max_y, max_z);
}


template<class T>
tBox3<T> tBox3<T>::make_union(const tBox3 &b1, const tBox3 &b2)
{
	return tBox3(min(b1.x_min(), b2.x_min()),
		    min(b1.y_min(), b2.y_min()),
		    min(b1.z_min(), b2.z_min()),
		    max(b1.x_max(), b2.x_max()),
		    max(b1.y_max(), b2.y_max()),
		    max(b1.z_max(), b2.z_max()));
}

template <class T>
inline void interval_intersection
(T mn_1, T mx_1, T mn_2, T mx_2, T &mn, T &mx)
{
	mn = std::max(mn_1, mn_2);
	mx = std::min(mx_1, mx_2);
}

template <class T>
tBox3<T> tBox3<T>::make_intersection(const tBox3 &b1, const tBox3 &b2)
{
	tBox3<T> result;
	for (size_t i = 0; i < 3; ++i) {
		T mn, mx;
		interval_intersection(b1.min_point()[i],
				      b1.max_point()[i],
				      b2.min_point()[i],
				      b2.max_point()[i], mn, mx);
		result.set_min(i, mn);
		result.set_min(i, mx);
	}
	return result;
}


/* Copyright (c) Hewlett Packard Company, 1997.  */

/* This program is freely distributable without licensing fees 
   and is provided without guarantee or warrantee expressed or 
   implied. This program is -not- in the public domain. */

#define GL_OCCLUSION_TEST_HP			0x8165
#define GL_OCCLUSION_RESULT_HP			0x8166
#define GL_OCCLUSION_TEST_HP_OLD		0x816E
#define GL_OCCLUSION_RESULT_HP_OLD		0x816F


#if 0
template<class T>
bool tBox3<T>::is_visible() const
{
#ifndef HP_EXT
	return true;
#else
	GLfloat x1 = _min_pt[0],
		y1 = _min_pt[1],
		z1 = _min_pt[2],
		x2 = _max_pt[0],
		y2 = _max_pt[1],
		z2 = _max_pt[2];
	GLboolean   result;

	glDisable(GL_STENCIL_TEST);
	glEnable(GL_OCCLUSION_TEST_HP);
	glDepthMask(GL_FALSE);
	glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);

	glBegin(GL_QUADS);

	/* Front */
	glVertex3f(x1,y1,z1);
	glVertex3f(x1,y2,z1);
	glVertex3f(x2,y2,z1);
	glVertex3f(x2,y1,z1);

	/* Right */
	glVertex3f(x2,y2,z1);
	glVertex3f(x2,y1,z1);
	glVertex3f(x2,y1,z2);
	glVertex3f(x2,y2,z2);

	/* Back */
	glVertex3f(x2,y2,z2);
	glVertex3f(x2,y1,z2);
	glVertex3f(x1,y1,z2);
	glVertex3f(x1,y2,z2);

	/* Left */
	glVertex3f(x1,y2,z2);
	glVertex3f(x1,y1,z2);
	glVertex3f(x1,y1,z1);
	glVertex3f(x1,y2,z1);

	/* Top */
	glVertex3f(x2,y2,z1);
	glVertex3f(x2,y2,z2);
	glVertex3f(x1,y2,z2);
	glVertex3f(x1,y2,z1);

	/* Bottom */
	glVertex3f(x2,y1,z1);
	glVertex3f(x2,y1,z2);
	glVertex3f(x1,y1,z2);
	glVertex3f(x1,y1,z1);
	glEnd();

	glGetBooleanv(GL_OCCLUSION_RESULT_HP, &result);
	glDisable(GL_OCCLUSION_TEST_HP);
	glEnable(GL_LIGHTING);
	glDepthMask(GL_TRUE);
	glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
	glEnable(GL_STENCIL_TEST);
	fprintf(stderr, "visible: %d\n", result);
	return (result == GL_TRUE);
#endif
}
#endif

/*
 * makes the box a little larger.
 * s - percent of the box width to add (can be negative > -0.5)
 */
template<class T>
void tBox3<T>::enlarge(value_type s)
{
    assert( (s >= 0) && (s <= 1));
    for (int i = 0; i < 3 ; ++i)
    {
        value_type w = _max_pt[i] - _min_pt[i];
        value_type addition = w * s;
        _min_pt[i] -= addition;
        _max_pt[i] += addition;
    }
}


template class tBox3<float>;
template class tBox3<double>;

GTB_END_NAMESPACE
