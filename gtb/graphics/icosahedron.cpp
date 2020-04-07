
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
#include <gtb/graphics/icosahedron.hpp>
#include <gtb/math/math.hpp>
#include <gtb/error/error.hpp>
#endif // WIN32

#ifdef OUTLINE
#define inline
#include <gtb/graphics/icosahedron.ipp>
#undef inline
#endif

using namespace std;

GTB_BEGIN_NAMESPACE


const double icosahedron::m_X = 0.525731112119133606;
const double icosahedron::m_Z = 0.850650808352039932;


GLfloat icosahedron::m_vdata[NUM_VERTICES][3] = {    
	{-m_X, 0.0, m_Z}, {m_X, 0.0, m_Z},
	{-m_X, 0.0, -m_Z}, {m_X, 0.0, -m_Z},    
	{0.0, m_Z, m_X}, {0.0, m_Z, -m_X},
	{0.0, -m_Z, m_X}, {0.0, -m_Z, -m_X},    
	{m_Z, m_X, 0.0}, {-m_Z, m_X, 0.0},
	{m_Z, -m_X, 0.0}, {-m_Z, -m_X, 0.0} 
};

GLuint icosahedron::m_tindices[NUM_FACES][3] = { 
	{0,4,1}, {0,9,4}, {9,5,4}, {4,5,8}, {4,8,1},    
	{8,10,1}, {8,3,10}, {5,3,8}, {5,2,3}, {2,7,3},    
	{7,10,3}, {7,6,10}, {7,11,6}, {11,0,6}, {0,1,6}, 
	{6,1,10}, {9,0,11}, {9,11,2}, {9,2,5}, {7,2,11}
};


icosahedron::icosahedron()
{
}


static void draw_triangle(GLfloat v1[3], GLfloat v2[3], GLfloat v3[3])
{
	glBegin(GL_TRIANGLES);
	glNormal3fv(v1);
	glVertex3fv(v1);
	glNormal3fv(v2);
	glVertex3fv(v2);
	glNormal3fv(v3);
	glVertex3fv(v3);
	glEnd();
}


static void normalize(GLfloat v[3])
{
	GLfloat d = sqrt((v[0] * v[0]) + (v[1] * v[1]) + (v[2] * v[2]));
	if (d == 0.0) {
		GTB_ERROR("zero length vector");
	}
	v[0] /= d;
	v[1] /= d;
	v[2] /= d;
}


static void render_subdivide(GLfloat v1[3],
			     GLfloat v2[3],
			     GLfloat v3[3],
			     int depth)
{
	if (depth == 0) {
		draw_triangle(v1, v2, v3);
		return;
	}
	GLfloat v12[3], v23[3], v31[3];
	for (int i = 0; i < 3; i++) {
		v12[i] = (v1[i] + v2[i]) / 2.0;
		v23[i] = (v2[i] + v3[i]) / 2.0;
		v31[i] = (v3[i] + v1[i]) / 2.0;
	}
	normalize(v12);
	normalize(v23);
	normalize(v31);
	render_subdivide(v1, v12, v31, depth - 1);
	render_subdivide(v2, v23, v12, depth - 1);
	render_subdivide(v3, v31, v23, depth - 1);
	render_subdivide(v12, v23, v31, depth - 1);
}


void icosahedron::render(unsigned depth)
{
	for (unsigned i = 0; i < NUM_FACES; i++) {    
		render_subdivide(m_vdata[m_tindices[i][0]],
				 m_vdata[m_tindices[i][1]],
				 m_vdata[m_tindices[i][2]],
				 depth);
	}
}


static void sample_points_subdivide(GLfloat v1[3],
				    GLfloat v2[3],
				    GLfloat v3[3],
				    vector<Point3> &points,
				    int depth)
{
	if (depth == 0) {
		return;
	}
	GLfloat v12[3], v23[3], v31[3];
	for (int i = 0; i < 3; i++) {
		v12[i] = (v1[i] + v2[i]) / 2.0;
		v23[i] = (v2[i] + v3[i]) / 2.0;
		v31[i] = (v3[i] + v1[i]) / 2.0;
	}
	normalize(v12);
	normalize(v23);
	normalize(v31);
	points.push_back(Point3(v12[0], v12[1], v12[2]));
	points.push_back(Point3(v23[0], v23[1], v23[2]));
	points.push_back(Point3(v31[0], v31[1], v31[2]));
	sample_points_subdivide(v1, v12, v31, points, depth - 1);
	sample_points_subdivide(v2, v23, v12, points, depth - 1);
	sample_points_subdivide(v3, v31, v23, points, depth - 1);
	sample_points_subdivide(v12, v23, v31, points, depth - 1);
}


void icosahedron::sample_points(std::vector<Point3> &points,
				unsigned depth)
{
	assert(points.size() == 0);
	for (unsigned i = 0; i < NUM_VERTICES; i++) {
		points.push_back(Point3(m_vdata[i][0],
					m_vdata[i][1],
					m_vdata[i][2]));
	}
	for (unsigned i = 0; i < NUM_FACES; i++) {    
		sample_points_subdivide(m_vdata[m_tindices[i][0]],
					m_vdata[m_tindices[i][1]],
					m_vdata[m_tindices[i][2]],
					points,
					depth);
	}
}


GTB_END_NAMESPACE
