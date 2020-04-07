
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
#include "extinction_coef_set.hpp"
#include <gtb/graphics/icosahedron.hpp>
#include <gtb/io/io.hpp>
#include <sstream>

#ifdef OUTLINE
#define inline
#include "extinction_coef_set.ipp"
#undef inline
#endif

using namespace std;

GTB_BEGIN_NAMESPACE


vector<Point3> extinction_coef_set::m_points;
bool extinction_coef_set::m_points_initialized =
	extinction_coef_set::init_points();


bool extinction_coef_set::init_points()
{
	icosahedron ico;
	unsigned depth = 0;	// 0: 12 points, 1: 72 points, 2: 312 points
	ico.sample_points(m_points, depth);
	return true;
}


extinction_coef_set::extinction_coef_set()
{
}


extinction_coef_set::~extinction_coef_set()
{
}


void extinction_coef_set::init(const OctreeNode &node)
{
	assert(m_points.size() > 0);
	for (unsigned i = 0; i < m_points.size(); i++) {
		Vector3 v = Vector3(m_points[i]);
		m_samples.push_back(extinction_coef(node, v));
	}
}


string extinction_coef_set::get_file_name(const char *model_file_name, int id)
{
	ostringstream ost;
	ost << model_file_name << "." << id << ".k";
	return ost.str();
}


void extinction_coef_set::read(const char *model_file_name, int id)
{
	assert(m_samples.size() == 0);
	string s = get_file_name(model_file_name, id);
	FILE *fp = xfopen(s.c_str(), "rb");
	unsigned n;
	read_unsigned(&n, fp);
	for (unsigned i = 0; i < n; i++) {
		float f;
		read_float(&f, fp);
		m_samples.push_back(extinction_coef(f));
	}
	assert(m_samples.size() == n);
	fclose(fp);
}


void extinction_coef_set::write(const char *model_file_name, int id)
{
	assert(m_samples.size() > 0);
	string s = get_file_name(model_file_name, id);
	FILE *fp = xfopen(s.c_str(), "wb");
	write_unsigned(m_samples.size(), fp);
	for (unsigned i = 0; i < m_samples.size(); i++) {
		write_float(m_samples[i], fp);
	}
	fclose(fp);
}


extinction_coef extinction_coef_set::closest(const Vector3 &v) const
{
	assert(m_points.size() >= 3);
	assert(m_samples.size() == m_points.size());

	// compute dot products
	vector<float> d(m_points.size());
	for (unsigned i = 0; i < d.size(); i++) {
		d[i] = v.dot(m_points[i]);
	}

	// find the three largest ones
	unsigned a, b, c;
	if (d[1] > d[0]) {
		a = 1;
		b = 0;
	} else {
		a = 0;
		b = 1;
	}
	if (d[2] > d[b]) {
		if (d[2] > d[a]) {
			c = b;
			b = a;
			a = 2;
		} else {
			c = b;
			b = 2;
		}
	} else {
		c = 2;
	}

	assert(a != b);
	assert(a != c);
	assert(b != c);
	assert(d[a] >= d[b]);
	assert(d[b] >= d[c]);

	for (unsigned i = 3; i < m_points.size(); i++) {
		if (d[i] > d[c]) {
			if (d[i] > d[b]) {
				if (d[i] > d[a]) {
					c = b;
					b = a;
					a = i;
				} else {
					c = b;
					b = i;
				}
			} else {
				c = i;
			}
		}
	}

	assert(a != b);
	assert(a != c);
	assert(b != c);
	assert(d[a] >= d[b]);
	assert(d[b] >= d[c]);
	for (unsigned i = 0; i < m_points.size(); i++) {
		if ((i == a) || (i == b) || (i == c)) {
			continue;
		}
		assert(d[c] >= d[i]);
	}

	// find weights
	Triangle3 t(m_points[a], m_points[b], m_points[c]);
	Point3 p(v.x(), v.y(), v.z());
	real_t wa, wb, wc;
	t.get_barycentric_coordinates(p, wa, wb, wc);

	// interpolate
	float f = ((wa * m_samples[a]) +
		   (wb * m_samples[b]) +
		   (wc * m_samples[c]));
	return extinction_coef(f);
}


GTB_END_NAMESPACE
