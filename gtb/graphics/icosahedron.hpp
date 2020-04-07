
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


#ifndef GTB_ICOSAHEDRON_INCLUDED
#define GTB_ICOSAHEDRON_INCLUDED

#include <gtb/common.hpp>
#include <gtb/graphics/box3.hpp>
#include <gtb/graphics/ogltools.h>

GTB_BEGIN_NAMESPACE


class icosahedron {
public:
	icosahedron();
	static void render(unsigned depth = 0);
	static Box3 bounding_box();
	static Point3 centroid();

	static void sample_points(std::vector<Point3> &points,
				  unsigned depth = 0);

	enum {
		NUM_VERTICES = 12,
		NUM_FACES = 20
	};

protected:
	static const double m_X;
	static const double m_Z;
	static GLfloat m_vdata[NUM_VERTICES][3];
	static GLuint m_tindices[NUM_FACES][3];
};


GTB_END_NAMESPACE

#ifndef OUTLINE
#include <gtb/graphics/icosahedron.ipp>
#endif

#endif // GTB_ICOSAHEDRON_INCLUDED
