
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

#include <gtb/graphics/box2.hpp>
#include <gtb/graphics/line2.hpp>
#include <gtb/graphics/segment2.hpp>
#include <gtb/graphics/ogltools.h>




#ifdef OUTLINE
#define inline
#include <gtb/graphics/box2.ipp>
#undef inline
#endif


using namespace std;


GTB_BEGIN_NAMESPACE


// reverse engineered from box3::_vertex_indices
template<class T>
const unsigned tBox2<T>::_vertex_indices[4][2] = {
	{0, 1},
	{1, 3},
	{3, 2},
	{2, 0}
};



template<class T>
typename tBox2<T>::value_type tBox2<T>::shortest_axis_length() const
{
	value_type dx = x_length();
	value_type dy = y_length();
	return min(dx, dy);
}


template<class T>
typename tBox2<T>::value_type tBox2<T>::longest_axis_length() const
{
	value_type dx = x_length();
	value_type dy = y_length();
	return max(dx, dy);
}


template<class T>
int tBox2<T>::classify_position(const tBox2 &b) const
{
	if ((_max_pt[0] < b._min_pt[0]) || (_min_pt[0] > b._max_pt[0]) ||
	    (_max_pt[1] < b._min_pt[1]) || (_min_pt[1] > b._max_pt[1])) {
		return OUTSIDE;
	}

	if ((_min_pt[0] <= b._min_pt[0]) && (_max_pt[0] >= b._max_pt[0]) &&
	    (_min_pt[1] <= b._min_pt[1]) && (_max_pt[1] >= b._max_pt[1])) {
		return INSIDE;
        }

	return INTERSECT;
}


template<class T>
tSegment2<T> tBox2<T>::face(unsigned i) const
{
	assert(i < 4);
	return tSegment2<T>(vertex(i, 0),
			    vertex(i, 1));
}


template<class T>
tLine2<T> tBox2<T>::plane(unsigned i) const
{
	assert(i < 4);
	return tLine2<T>(vertex(_vertex_indices[i][0]),
		     vertex(_vertex_indices[i][1]));
}


template<class T>
void tBox2<T>::render() const
{
	glBegin(GL_QUADS);
	for (unsigned i = 0; i < 4; i++) {
		vertex(i, 0).load();
	}
	glEnd();
}


template<class T>
void tBox2<T>::outline() const
{
	glPushAttrib(GL_ENABLE_BIT);
	glDisable(GL_LIGHTING);

	glBegin(GL_LINE_LOOP);
	vertex(0).load();
	vertex(1).load();
	vertex(2).load();
	vertex(3).load();
	glEnd();

	glPopAttrib();
}

template<class T>
tBox2<T> tBox2<T>::bounding_box(const vector<tPoint2<T> > &v)
{
	assert(v.size() > 0);
	tPoint2<T> p = v[0];
	value_type min_x = p.x(), min_y = p.y();
	value_type max_x = p.x(), max_y = p.y();
	for (unsigned i = 1; i < v.size(); i++) {
		p = v[i];
		min_x = min(min_x, p.x());
		min_y = min(min_y, p.y());
		max_x = max(max_x, p.x());
		max_y = max(max_y, p.y());
	}
	return tBox2(min_x, min_y,
 		     max_x, max_y);
}


template<class T>
tBox2<T> tBox2<T>::make_union(const tBox2 &b1, const tBox2 &b2)
{
	return tBox2(min(b1.x_min(), b2.x_min()),
		    min(b1.y_min(), b2.y_min()),
		    max(b1.x_max(), b2.x_max()),
		    max(b1.y_max(), b2.y_max()));
}


/* Copyright (c) Hewlett Packard Company, 1997.  */

/* This program is freely distributable without licensing fees 
   and is provided without guarantee or warrantee expressed or 
   implied. This program is -not- in the public domain. */

#define GL_OCCLUSION_TEST_HP			0x8165
#define GL_OCCLUSION_RESULT_HP			0x8166
#define GL_OCCLUSION_TEST_HP_OLD		0x816E
#define GL_OCCLUSION_RESULT_HP_OLD		0x816F


/*
 * makes the box a little larger.
 * s - percent of the box width to add (can be negative > -0.5)
 */
template<class T>
void tBox2<T>::enlarge(value_type s)
{
    assert( (s >= 0) && (s <= 1));
    for (int i = 0; i < 2 ; ++i)
    {
        value_type w = _max_pt[i] - _min_pt[i];
        value_type addition = w * s;
        _min_pt[i] -= addition;
        _max_pt[i] += addition;
    }
}

GTB_GENERATE_CLASS_TYPEDEFS(Box2)

template class tBox2<float>;
template class tBox2<double>;

GTB_END_NAMESPACE
