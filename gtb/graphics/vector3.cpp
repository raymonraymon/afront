
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
#include <gtb/graphics/vector3.hpp>
#include <gtb/graphics/point3.hpp>
#include <gtb/math/math.hpp>
#endif // WIN32


#ifdef OUTLINE
#define inline
#include <gtb/graphics/vector3.ipp>
#undef inline
#endif


GTB_BEGIN_NAMESPACE

template<>
const tVector3<float> tVector3<float>::VECTOR3_ZERO(0.0, 0.0, 0.0);
template<>
const tVector3<float> tVector3<float>::VECTOR3_POSITIVE_X(1.0, 0.0, 0.0);
template<>
const tVector3<float> tVector3<float>::VECTOR3_NEGATIVE_X(-1.0, 0.0, 0.0);
template<>
const tVector3<float> tVector3<float>::VECTOR3_POSITIVE_Y(0.0, 1.0, 0.0);
template<>
const tVector3<float> tVector3<float>::VECTOR3_NEGATIVE_Y(0.0, -1.0, 0.0);
template<>
const tVector3<float> tVector3<float>::VECTOR3_POSITIVE_Z(0.0, 0.0, 1.0);
template<>
const tVector3<float> tVector3<float>::VECTOR3_NEGATIVE_Z(0.0, 0.0, -1.0);

template<>
const tVector3<double> tVector3<double>::VECTOR3_ZERO(0.0, 0.0, 0.0);
template<>
const tVector3<double>  tVector3<double>::VECTOR3_POSITIVE_X(1.0, 0.0, 0.0);
template<>
const tVector3<double>  tVector3<double>::VECTOR3_NEGATIVE_X(-1.0, 0.0, 0.0);
template<>
const tVector3<double>  tVector3<double>::VECTOR3_POSITIVE_Y(0.0, 1.0, 0.0);
template<>
const tVector3<double>  tVector3<double>::VECTOR3_NEGATIVE_Y(0.0, -1.0, 0.0);
template<>
const tVector3<double>  tVector3<double>::VECTOR3_POSITIVE_Z(0.0, 0.0, 1.0);
template<>
const tVector3<double>  tVector3<double>::VECTOR3_NEGATIVE_Z(0.0, 0.0, -1.0);

template<>
tVector3<double>& tVector3<double>::rotate(const tVector3 &axis, value_type theta)
{
	// From Goldstein
	*this =	*this * cos(theta) +
		axis * this->dot(axis) * (1.0 - cos(theta)) -
		this->cross(axis) * sin(theta);
	return *this;
}

template<>
tVector3<float>& tVector3<float>::rotate(const tVector3 &axis, value_type theta)
{
	// From Goldstein
	*this =	*this * cosf(theta) +
		axis * this->dot(axis) * (1.0f - cosf(theta)) -
		this->cross(axis) * sinf(theta);
	return *this;
}

// Finds the rotation that takes V to W (Graphics Gems I, p. 531).
// (Here we pre-multiply, though, i.e., M * p, so we have to use the
// transposetemplate<class T> of what's shown there.)
template<class T>
tMatrix4<T> tVector3<T>::rotation(const tVector3 &V, const tVector3 &W)
{
	tVector3 N = V.cross(W);
	N.normalize();

	tVector3 M = N.cross(V);
	M.normalize();

	tMatrix4<T> Q(V.x(), V.y(), V.z(), 0.0,
		  M.x(), M.y(), M.z(), 0.0,
		  N.x(), N.y(), N.z(), 0.0,
		  0.0, 0.0, 0.0, 1.0);

	tVector3 W2 = Q * W;

	tMatrix4<T> RT(W2.x(), -W2.y(), 0.0, 0.0,
		   W2.y(), W2.x(), 0.0, 0.0,
		   W2.z(), W2.z(), 1.0, 0.0,
		   0.0, 0.0, 0.0, 1.0);

	return Q.transpose() * RT * Q;
}

template<class T>
void tVector3<T>::element_wise_multiply(const tVector3& rhs)
{
    _v[0] *= rhs[0];
    _v[1] *= rhs[1];
    _v[2] *= rhs[2];
}


template class tVector3<float>;
template class tVector3<double>;


GTB_END_NAMESPACE
