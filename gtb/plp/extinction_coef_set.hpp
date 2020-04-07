
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


#ifndef EXTINCTION_COEF_SET_INCLUDED
#define EXTINCTION_COEF_SET_INCLUDED

#include <vector>
#include <string>
#include "extinction_coef.hpp"

GTB_BEGIN_NAMESPACE


class extinction_coef_set {
public:
	extinction_coef_set();
	~extinction_coef_set();
	void init(const OctreeNode &node);
	void read(const char *model_file_name, int id);
	void write(const char *model_file_name, int id);
	extinction_coef closest(const Vector3 &v) const;

protected:
	static std::string get_file_name(const char *model_file_name, int id);
	static bool init_points();

	std::vector<extinction_coef> m_samples;
	static std::vector<Point3> m_points;
	static bool m_points_initialized;
};


GTB_END_NAMESPACE

#ifndef OUTLINE
#include "extinction_coef_set.ipp"
#endif

#endif // EXTINCTION_COEF_SET_INCLUDED
