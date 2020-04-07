
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


#ifndef __RSTPLANE
#define __RSTPLANE

#include <gtb/common.hpp>
#include <gtb/real/real.hpp>
#include <gtb/math/math.hpp>
#include <gtb/graphics/vector3.hpp>

GTB_BEGIN_NAMESPACE

/*
 * This class represents a plane in 3D
 * using two angles and a distance from the origin
 */
template <class T>
struct tCPlaneRST
{
    typedef T value_type;

    tCPlaneRST() {}
    tCPlaneRST(value_type v_r, value_type v_s, value_type v_t) : r(v_r), s(v_s), t(v_t) {}

    value_type r,s; // Angles
    value_type t;   // Distance from the origin
};

typedef tCPlaneRST<float> CPlaneRSTf;
typedef tCPlaneRST<double> CPlaneRSTd;
#if defined(REAL_IS_FLOAT)
typedef CPlaneRSTf CPlaneRST;
#else
typedef CPlaneRSTd CPlaneRST;
#endif


GTB_END_NAMESPACE

#endif // __RSTPLANE
