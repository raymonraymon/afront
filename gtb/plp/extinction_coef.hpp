
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


#ifndef EXTINCTION_COEF_INCLUDED
#define EXTINCTION_COEF_INCLUDED

#include <gtb/graphics/octree_node.hpp>

GTB_BEGIN_NAMESPACE


class extinction_coef {
public:
	explicit extinction_coef(float value);
	explicit extinction_coef(const OctreeNode &node);
	extinction_coef(const OctreeNode &node, const Vector3 &v);
	operator float() const;

protected:
	float m_value;
};


GTB_END_NAMESPACE

#ifndef OUTLINE
#include "extinction_coef.ipp"
#endif

#endif // EXTINCTION_COEF_INCLUDED
